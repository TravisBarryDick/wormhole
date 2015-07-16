/**
 * @file  rbcBinary _parser.h
 * @brief  parse row block containers in binary format
 */
#pragma once
#include <limits>
#include "data/row_block.h"
#include "data/parser.h"

namespace dmlc {
namespace data {

/**
 * \brief RowBlockContainer Binary Parser
 */
template <typename IndexType>
class RowBlockContainerBinaryParser : public ParserImpl<IndexType> {
 public:
  explicit RowBlockContainerBinaryParser(const char *filename)
      : bytes_read_(0), filename_(filename) { 
    //open file
    outfile_ = Stream::Create(filename, "r");
  }
  virtual ~RowBlockContainerBinaryParser () {
    delete outfile_;
  }

  virtual void BeforeFirst(void) {}
  virtual size_t BytesRead(void) const {
    return bytes_read_;
  }
  virtual bool ParseNext(std::vector<RowBlockContainer<IndexType> > *data) {
  
    
    data->resize(1);
    RowBlockContainer<IndexType>& blk = (*data)[0];
    blk.Clear();
    bool read = blk.Load(outfile_);
    if (!read) return false;
    bytes_read_ += blk.MemCostBytes();
    return true;
  }

 private:
  // number of bytes readed
  size_t bytes_read_;
  const char *filename_;
  Stream *outfile_;
};

}  // namespace data
}  // namespace dmlc
