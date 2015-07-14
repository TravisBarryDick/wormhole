#pragma once
#include "dmlc/data.h"
#include "data/row_block.h"
#include "dmlc/io.h"
#include <string>

namespace dddml{
using RowBlockContainer = dmlc::data::RowBlockContainer;

template <typename I>
public class BufferedWriter{
	
	private:
	dmlc::Stream *out_stream;
	std::string filename;
	RowBlockContainer<I> buffer;
	size_t max_buffer_size;
	bool initialized;
	void Flush();
	
	public:
	BufferedWriter(const char *filename, size_t buffer_size = 1000000) : filename(filename), max_buffer_size(buffer_size) , initialized(false)
	{
		buffer.Clear();
	}
		
	void Write(const Row<I> &row);
}
}//namespace dddml
