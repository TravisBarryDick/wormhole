#pragma once
#include "dmlc/data.h"
#include "data/row_block.h"
#include "dmlc/io.h"
#include <string>

namespace dddml{

template <typename I>
class BufferedWriter{

	private:
	dmlc::Stream *out_stream;
	std::string filename;
	dmlc::data::RowBlockContainer<I> buffer;
	size_t max_buffer_size;
	bool initialized;

	public:
	void Flush();

	BufferedWriter(const char *filename, size_t buffer_size = 1000000) : filename(filename), max_buffer_size(buffer_size) , initialized(false)
	{
		buffer.Clear();
	}

	void Write(const dmlc::Row<I> &row);
};

template <typename I>
inline void BufferedWriter<I>::Flush()
{
	char flags[2] = {'\0', '\0'};
	flags[0] = (initialized) ? 'a' : 'w' ;
	initialized = true;
	out_stream = dmlc::Stream::Create(filename.c_str(), flags);
	buffer.Save(out_stream);
	delete out_stream;
	buffer.Clear();
}

template <typename I>
inline void BufferedWriter<I>::Write(const dmlc::Row<I> &row)
{
	buffer.Push(row);
	if (buffer.Size() >= max_buffer_size) {
		Flush();
	}
}

}//namespace dddml
