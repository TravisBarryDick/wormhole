#pragma once
#include "dmlc/data.h"
#include "data/row_block.h"
#include "dmlc/io.h"
#include <string>

namespace dddml{

template <typename I>
class BufferedWriter{

	private:
	dmlc::ostream *out;
	std::string filename;

	public:
  void Flush(){}

	BufferedWriter(const char *filename) : filename(filename) 
  {
	  //char flags[2] = {'\0', '\0'};
	  //flags[0] = (initialized) ? 'a' : 'w' ;
	  auto out1 = dmlc::Stream::Create(filename, "w");
    out = new dmlc::ostream(out1);
  }

  ~BufferedWriter(){delete out;}

	void Write(const dmlc::Row<I> &row);
};


template <typename I>
inline void BufferedWriter<I>::Write(const dmlc::Row<I> &row)
{
    *out << row.label ; //label
    for (int j = 0; j < row.length; ++j)
    {
      *out << ' ' << row.index[j] ;
      if (row.value) *out  << ':' << row.value[j];
    }
    *out << '\n';
}

}//namespace dddml
