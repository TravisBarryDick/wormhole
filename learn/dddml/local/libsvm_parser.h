#pragma once 

#include <vector>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cerrno>
#include <cassert>
#include "data.h"
#include "./row_block.h"
#include "str2num.h"
//debug
// #include <iostream>
//#include <cstdio>

namespace dmlc{
namespace data{


template <typename IndexType>
inline void 
ParseBlock(char *begin,
           char *end,
           RowBlockContainer<IndexType> *out) {
  std::cerr << "******\nNOTE: Check indexing. Currently using indexing starting from 1.!!!!!!!\n***********\n";
  out->Clear();
  char * lbegin = begin;
  char * lend = lbegin;
  while (lbegin != end) {
    // get line end
    lend = lbegin + 1;
    while (lend != end && *lend != '\n' && *lend != '\r') ++lend;
    // parse label[:weight]
    const char * p = lbegin;
    const char * q = NULL;
    real_t label;
    real_t weight;
    int r = ParsePair<real_t, real_t>(p, lend, &q, label, weight);

    if (r < 1) {
      // empty line
      lbegin = lend;
      continue;
    }
    if (r == 2) {
      // has weight
      out->weight.push_back(weight);
    }
    if (out->label.size() != 0) {
      out->offset.push_back(out->index.size());
    }
    out->label.push_back(label);
    // parse feature[:value]
    p = q;
    while (p != lend) {
      IndexType featureId;
      real_t value;
      int r = ParsePair<IndexType, real_t>(p, lend, &q, featureId, value);
      if (r < 1) {
        p = q;
        continue;
      }
      --featureId;//!!!!!!IMPORTANT: For zero based indexing ///////////TODO///////////////!!!!!!!!!!!
      //std::cerr << "******\nNOTE: Check indexing. Currently using indexing starting from 1.!!!!!!!\n***********\n";
      out->index.push_back(featureId);
      if (r == 2) {
        // has value
        out->value.push_back(value);
      }
      p = q;
    }
    // next line
    lbegin = lend;
  }
  if (out->label.size() != 0) {
    out->offset.push_back(out->index.size());
  }
  assert(out->label.size() + 1 == out->offset.size());
}


inline std::string get_file_contents(const char *filename)
{
  std::ifstream in(filename, std::ios::in | std::ios::binary);
  if (!in.is_open())
    std::cout << "FILE NOT FOUND\n";
  if (in)
  {
    std::ostringstream contents;
    contents << in.rdbuf();
    in.close();
    return(contents.str());
  }
  throw(errno);
}

// trim from left
inline std::string& ltrim(std::string& s, const char* t = " \t\n\r\f\v")
{
    s.erase(0, s.find_first_not_of(t));
    return s;
}

// trim from right
inline std::string& rtrim(std::string& s, const char* t = " \t\n\r\f\v")
{
    s.erase(s.find_last_not_of(t) + 1);
    return s;
}

// trim from left & right
inline std::string& trim(std::string& s, const char* t = " \t\n\r\f\v")
{
    return ltrim(rtrim(s, t), t);
}

// copying versions

inline std::string ltrim_copy(std::string s, const char* t = " \t\n\r\f\v")
{
    return ltrim(s, t);
}

inline std::string rtrim_copy(std::string s, const char* t = " \t\n\r\f\v")
{
    return rtrim(s, t);
}

inline std::string trim_copy(std::string s, const char* t = " \t\n\r\f\v")
{
    return trim(s, t);
}


RowBlockContainer<int>  libsvmread(const char *filename)
{
  std::string contents = get_file_contents(filename);
  //std::cout << contents << std::endl;
  trim(contents);
  char * cstr = new char [contents.length()+1];
  std::strcpy (cstr, contents.c_str());
  char *start = cstr;
  char *finish = start + contents.length() - 1;

  RowBlockContainer<int> rbc;

  //std::cout << "done 1\n";

  ParseBlock(start, finish, &rbc);

  //std::cout << "done 2\n";

  return rbc;
  throw(errno);
}

template <typename I>
void libsvmwrite(const char *filename, RowBlockContainer<I> &rbc)
{
  int ptr = 0;
  std::ofstream out;
  out.open(filename);
  for (int i = 0; i < rbc.label.size(); ++i )
  {
    out << rbc.label[i]; //label
    for (int j = rbc.offset[i]; j < rbc.offset[i+1]; ++j)
    {
      out << ' ' << rbc.index[j] << ':' << rbc.value[j];
    }
    out << std::endl;
  }
  out.close();
}

template <typename I>
void libsvmwrite(const char *filename, RowBlock<I> &rbc)
{
  int ptr = 0;
  std::ofstream out;
  out.open(filename);
  for (int i = 0; i < rbc.size; ++i )
  {
    out << rbc.label[i]; //label
    for (int j = rbc.offset[i]; j < rbc.offset[i+1]; ++j)
    {
      out << ' ' << rbc.index[j] << ':' << rbc.value[j];
    }
    out << std::endl;
  }
  out.close();
}

template <typename I>
void libsvmwrite(const char *filename, const RowBlock<I> &rbc, std::vector<int> &assignments)
{
  int ptr = 0;
  std::ofstream out;
  out.open(filename);
  for (int i = 0; i < rbc.size; ++i )
  {
    out << assignments[i]; //label
    for (int j = rbc.offset[i]; j < rbc.offset[i+1]; ++j)
    {
      out << ' ' << rbc.index[j] << ':' << rbc.value[j];
    }
    out << std::endl;
  }
  out.close();
}

template <typename I>
void libsvmwrite(const char *filename, const RowBlock<I> &rbc, std::vector<std::vector<int>> &assignments)
{
  int ptr = 0;
  std::ofstream out;
  out.open(filename);
  assert(rbc.size == assignments.size());
  for (int i = 0; i < rbc.size; ++i )
  {
    //out << assignments[i]; //label
    for (int j = 0; j < assignments[i].size(); ++j)
    {
      out << assignments[i][j] << ((j == (assignments[i].size() - 1)) ? ' ' : ',' );
    }
    for (int j = rbc.offset[i]; j < rbc.offset[i+1]; ++j)
    {
      out << ' ' << rbc.index[j] << ':' << rbc.value[j];
    }
    out << std::endl;
  }
  out.close();
}

}
} //namespace


// int main()
// {
//   dmlc::data::RowBlockContainer<int> rbc = dmlc::data::libsvmread("./sample.txt");
//   std::cout << "REad \n";
//   dmlc::data::libsvmwrite("./sample2.txt", rbc);
// }






