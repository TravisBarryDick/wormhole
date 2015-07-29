/**
 * @file   async_sgd.h
 * @brief  A template to implement async sgd in parameter server
 */
#include "ps.h"
#include "ps/app.h"
#include "base/minibatch_iter.h"
#include "base/progress.h"
#include "base/dist_monitor.h"
#include "base/workload_pool.h"
#include "base/string_stream.h"

namespace dmlc {
namespace solver {

// commands
static const int kProcess = 1;
static const int kSaveModel = 2;
static const int kLoadModel = 3;
static const int kMaxNumCmd = 10;

/**************************************************************************
 * \brief The scheduler node, which issues workloads to workers/servers, in
 * charge of fault tolerance, and also prints the progress
 **************************************************************************/
template <typename Progress>
class AsyncSGDScheduler : public ps::App {
 protected:
  std::string train_data_;
  std::string val_data_;
  std::string data_format_;
  std::string model_in_;
  std::string model_out_;

  bool worker_local_data_ = false;
  int load_model_ = 0;
  int save_model_ = 0;
  int num_part_per_file_ = 10;
  int max_data_pass_ = 1;
  int cur_data_pass_ = 1;  // start from 1
  Workload::Type cur_type_;
  int disp_itv_ = 1;
  bool predict_ = false;
  virtual bool Stop(const Progress& cur, const Progress& prev, bool train) {
    return false;
  }

  template <typename Config>
  void Init(const Config& conf) {
    worker_local_data_ = conf.use_worker_local_data();
    train_data_        = conf.train_data();
    val_data_          = conf.val_data();
    data_format_       = conf.data_format();
    num_part_per_file_ = conf.num_parts_per_file();
    max_data_pass_     = conf.max_data_pass();
    disp_itv_          = conf.disp_itv();
    save_model_        = conf.save_model();
    load_model_        = conf.load_model();
    model_in_          = conf.model_in();
    model_out_         = conf.model_out();
    predict_           = conf.pred_out().size() > 0;
  }

 public:
  AsyncSGDScheduler() {
    sys_.manager().AddNodeFailureHandler([this](const std::string& id) {
        pool_.Reset(id);
      });
  }
  virtual ~AsyncSGDScheduler() { }

  virtual void ProcessResponse(ps::Message* response) {
    if (response->task.cmd() == kProcess) {
      auto id = response->sender;
      pool_.Finish(id);

      if (response->task.msg().size()) {
        CHECK(worker_local_data_);
        StringStream ss(response->task.msg());
        Workload wl; wl.Load(&ss);
        pool_.Add(wl.file, num_part_per_file_, id);
      }

      Workload wl; pool_.Get(id, &wl);
      if (!wl.Empty()) {
        CHECK_EQ(wl.file.size(), (size_t)1);
        wl.type = cur_type_;
        wl.data_pass = cur_data_pass_;
        wl.file[0].format = data_format_;
        SendWorkload(id, wl);
      }
    }
  }

  virtual bool Run() {
    printf("connected %d servers and %d workers\n",
           ps::NumServers(), ps::NumWorkers());
    start_time_ = GetTime();

    if (predict_) {
      CHECK(model_in_.size()) << "should provide model_in for prediction";
    }

    if (model_in_.size()) {
      if (load_model_ > 0) {
        printf("loading model from #iter = %d\n", load_model_);
        cur_data_pass_ = load_model_;
      } else {
        printf("loading the last model\n");
        cur_data_pass_ = 0;
      }
      ps::Task task; task.set_cmd(kLoadModel + load_model_ * kMaxNumCmd);
      Wait(Submit(task, ps::kServerGroup));
      Iterate(Workload::VAL);
      ++ cur_data_pass_;
    }

    if (predict_) {
      printf("prediction is done!\n");
      return true;
    }

    for (; cur_data_pass_ <= max_data_pass_; ++cur_data_pass_) {
      if (Iterate(Workload::TRAIN) || Iterate(Workload::VAL)) {
        printf("hit stop critera\n"); break;
      }
      if (cur_data_pass_ == max_data_pass_) {
        printf("hit max number of data passes\n");
      }
      SaveModel(false);
    }

    SaveModel(true);
    printf("training is done!\n");
    return true;
  }

 private:
  void SaveModel(bool force) {
    if (model_out_.size() == 0) return;
    if (force || (save_model_ > 0 && cur_data_pass_ % save_model_ == 0)) {
      int iter = force ? 0 : cur_data_pass_;
      if (iter == 0) {
        printf("saving final model to %s\n", model_out_.c_str());
      } else {
        printf("saving model to %s-iter_%d\n", model_out_.c_str(), iter);
      }
      ps::Task task; task.set_cmd(kSaveModel + iter * kMaxNumCmd);
      Wait(Submit(task, ps::kServerGroup));
    }
  }

  // return true if time for stop
  bool Iterate(Workload::Type type) {
    cur_type_ = type;
    bool stop = false;
    std::string data;
    bool is_train = type == Workload::TRAIN;
    if (is_train) {
      printf("training #iter = %d\n", cur_data_pass_);
      data = train_data_;
    } else {
      if (predict_) {
        CHECK(val_data_.size()) << "should provide val_data for prediction";
        printf("predicting\n");
      } else {
        printf("validating #iter = %d\n", cur_data_pass_);
      }
      data = val_data_;
    }
    if (data.empty()) return stop;

    pool_.Clear();
    if (!worker_local_data_) {
      Workload wl; pool_.Match(data, &wl);
      pool_.Add(wl.file, num_part_per_file_);
      if (predict_) {
        CHECK_EQ(wl.file.size(), (size_t)1)
            << "use single file for prediction";
      }
      size_t nwl = wl.file.size() * num_part_per_file_;
      if (cur_data_pass_ == 1 && (nwl < ps::NumWorkers())) {
        fprintf(stderr, "WARNING: # of data file (%d) < # of workers (%d)\n",
                (int)nwl, (int)ps::NumWorkers());
        fprintf(stderr, "         You may want to increase \"num_parts_per_file\"\n");
      }
    }

    // send an empty workerload to all workers
    Workload wl; wl.type = type; SendWorkload(ps::kWorkerGroup, wl);

    // print every k sec for training, while print at the end for validation
    printf("  sec %s\n", prog_.HeadStr().c_str());
    fflush(stdout);

    while (!pool_.IsFinished()) {
      sleep(disp_itv_);
      if (is_train) {
        if (ShowProgress(is_train)) {
          stop = true;
          pool_.ClearRemain();
        }
      }
    }
    if (type != Workload::TRAIN) {
      stop = ShowProgress(is_train);
    }
    return stop;
  }

  // return true if it's time for stopping
  bool ShowProgress(bool is_train) {
    bool ret = false;
    Progress cur;
    monitor_.Get(&cur); monitor_.Clear();
    auto disp = cur.PrintStr(&prog_);
    if (disp.empty()) return ret;
    printf("%5.0lf  %s\n", GetTime() - start_time_, disp.c_str());
    fflush(stdout);

    if (Stop(cur, prog_, is_train)) ret = true;
    prog_.Merge(&cur);
    return ret;
  }

  void SendWorkload(const std::string id, const Workload& wl) {
    StringStream ss; wl.Save(&ss);
    ps::Task task; task.set_msg(ss.str());
    task.set_cmd(kProcess);
    Submit(task, id);
  }

  double start_time_;
  WorkloadPool pool_;
  bool done_ = false;
  Progress prog_;
  ProgressMonitor<Progress> monitor_;
};

/**************************************************************************
 * \brief A server node, which maintains a part of the model, and update the
 * model once received gradients from workers
 **************************************************************************/
class AsyncSGDServer : public ps::App {
 protected:
  /**
   * \brief Report the progress to the scheduler
   */
  void Report(const IProgress* const prog) {
    reporter_.Report(prog);
  }

  virtual void SaveModel(int iter) = 0;
  virtual void LoadModel(int iter) = 0;

  std::string ModelName(const std::string& base, int iter) {
    CHECK(base.size()) << "empty model name";
    std::string name = base;
    if (iter > 0) name += "_iter-" + std::to_string(iter);
    return name + "_part-" + std::to_string(ps::MyRank());
  }

 public:
  AsyncSGDServer() {}
  virtual ~AsyncSGDServer() { }

  virtual void ProcessRequest(ps::Message* request) {
    int c = request->task.cmd();
    int cmd = c % kMaxNumCmd;
    int iter = c / kMaxNumCmd;

    if (cmd == kSaveModel) {
      SaveModel(iter);
    } else if (cmd == kLoadModel) {
      LoadModel(iter);
    }
  }

 private:
  ProgressReporter reporter_;

};

/*****************************************************************************
 * \brief A worker node, which takes a part of training data and calculate the
 * gradients
 *****************************************************************************/
class AsyncSGDWorker : public ps::App {
  /// for applications
 protected:
  // feature id type
  using FeaID = ps::Key;
  // a minibatch data
  using Minibatch = dmlc::RowBlock<FeaID>;

  /**
   * \brief Process one minibatch
   */
  virtual void ProcessMinibatch(const Minibatch& mb, const Workload& wl) = 0;

  /**
   * \brief Mark one minibatch is finished
   * one must call this function when the minibatch is actually done
   */
  void FinishMinibatch();

  /**
   * \brief Report the progress to the scheduler
   */
  void Report(const IProgress *const prog) {
    reporter_.Report(prog);
  }

  std::string PredictName(const std::string& base, int part) {
    CHECK(base.size());
    char* fname = new char[1000];
    snprintf(fname, 1000, "%s_part-%02d", base.c_str(), part);
    return std::string(fname);
  }

  int minibatch_size_ = 10000;
  int shuffle_ = 0;
  int max_delay_ = 0;

  // for validation test
  int val_minibatch_size_ = 1000000;
  int val_max_delay_ = 10;

  bool predict_ = false;

  bool worker_local_data_ = false;
  std::string train_data_;
  std::string val_data_;

  double workload_time_;

  float neg_sampling_ = 1.0;

  template <typename Config>
  void Init(const Config& conf) {
    minibatch_size_ = conf.minibatch();
    shuffle_        = conf.rand_shuffle();
    max_delay_      = conf.max_delay();
    predict_        = conf.pred_out().size() > 0;
    neg_sampling_   = conf.neg_sampling();
    if (conf.use_worker_local_data()) {
      train_data_        = conf.train_data();
      val_data_          = conf.val_data();
      worker_local_data_ = true;
    }
  }
  /// implement system APIs
 public:
  AsyncSGDWorker() { }
  virtual ~AsyncSGDWorker() { }

  // process request from the scheduler node
  virtual void ProcessRequest(ps::Message* request) {
    int cmd = request->task.cmd();
    if (cmd == kProcess) {
      StringStream ss(request->task.msg());
      Workload wl; wl.Load(&ss);
      if (wl.Empty()) {
        if (!worker_local_data_) return;

        Workload local;
        if (wl.type == Workload::TRAIN) {
          WorkloadPool::Match(train_data_, &local);
        } else {
          WorkloadPool::Match(val_data_, &local);
        }
        StringStream ss; local.Save(&ss);
        ps::Task res; res.set_msg(ss.str());
        Reply(request, res);
      } else {
        Process(wl);
      }
    }
  }

 private:
  void Process(const Workload& wl) {
    bool train = wl.type == Workload::TRAIN;
    int mb_size = train ? minibatch_size_ : val_minibatch_size_;
    int shuffle = train ? minibatch_size_ * shuffle_ : 0;
    int max_delay = predict_ ? 0 : (train ? max_delay_ : val_max_delay_);
    float neg_sampling = train ? neg_sampling_ : 1.0;
    LOG(INFO) << ps::MyNodeID() << ": " << wl.ShortDebugString()
              << ", minibatch = " << mb_size
              << ", max_delay = " <<  max_delay
              << ", shuffle size = " << shuffle
              << ", negative sampling = " << neg_sampling;

    num_mb_fly_ = num_mb_done_ = 0;
    start_ = GetTime();
    workload_time_ = 0;

    CHECK_EQ(wl.file.size(), (size_t)1);
    auto file = wl.file[0];
    dmlc::data::MinibatchIter<FeaID> reader(
        file.filename.c_str(), file.k, file.n, file.format.c_str(),
        mb_size, shuffle, neg_sampling);
    reader.BeforeFirst();
    while (reader.Next()) {
      // wait for data consistency
      WaitMinibatch(max_delay);

      ProcessMinibatch(reader.Value(), wl);

      mb_mu_.lock(); ++ num_mb_fly_; mb_mu_.unlock();
    }

    // wait untill all are done
    WaitMinibatch(0);
    LOG(INFO) << ps::MyNodeID() << ": finished";
  }

  // wait if the currenta number of on processing minibatch > num
  inline void WaitMinibatch(int num) {
    std::unique_lock<std::mutex> lk(mb_mu_);
    mb_cond_.wait(lk, [this, num] {return num_mb_fly_ <= num;});
  }

  int num_mb_fly_;
  int num_mb_done_;
  std::mutex mb_mu_;
  std::condition_variable mb_cond_;
  double start_;
  ProgressReporter reporter_;
};

void AsyncSGDWorker::FinishMinibatch() {
  // wake the main thread
  mb_mu_.lock();
  -- num_mb_fly_;
  ++ num_mb_done_;
  mb_mu_.unlock();
  mb_cond_.notify_one();

  // log some info

  double ttl_time = (GetTime() - start_);
  std::string overhead;
  if (workload_time_ > 0) {
    int avg_oh
        = std::max(ttl_time-workload_time_, (double)0) / ttl_time * 100;
    overhead = "overhead " + std::to_string(avg_oh) + "%, ";
  }
  LOG(INFO) << num_mb_done_ << " done, avg time "
            << ttl_time / num_mb_done_ << ", " << overhead
            << num_mb_fly_ << " on running";
}


}  // namespace solver
}  // namespace dmlc
