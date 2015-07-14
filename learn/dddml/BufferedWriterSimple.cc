#include "BufferedWriter.h"

namespace dddml{

template <typename I>
inline void BufferedWriter<I>::Flush()
{
	char flags[2] = {'\0', '\0'};
	flags[0] = (initialized) ? 'a' : 'w' ;
	intialized = true;
	out_stream = dmlc::Stream::Create(filename.c_str()), flags);	
	buffer.Save(out_stream);
	delete out_stream;
	buffer.Clear();
}

template <typename I>
inline void BufferedWriter<I>::Write(const Row<I> &row)
{
	buffer.push(row);
	if (buffer.Size() >= max_buffer_size) {
		Flush();
	}
}

} //namespace dddml
