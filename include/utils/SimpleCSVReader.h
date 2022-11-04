#pragma once

#include <string>
#include <vector>
#include <unordered_map>

// An improved CSV reader which reads in CSV files and doesn't cast values to strings. Required due
// to cross-section formatting in OpenMC.
class SimpleCSVReader
{
public:
  SimpleCSVReader(const std::string & file, const std::string & delimiter = ",");

  bool read();

  std::size_t numColumns() const { return _column_headers.size(); }

  std::string & getHeader(std::size_t column_index) { return _column_headers[column_index]; }

  std::vector<std::string> & getHeaders() { return _column_headers; }

  std::vector<std::string> & getColumn(const std::string & header) { return _csv_data[header]; }

  std::vector<std::string> & getColumn(std::size_t column_index)
  {
    return _csv_data[_column_headers[column_index]];
  }

  std::string & getEntry(const std::string & header, std::size_t row)
  {
    return _csv_data[header][row];
  }

  std::string & getEntry(std::size_t column_index, std::size_t row)
  {
    return _csv_data[_column_headers[column_index]][row];
  }

private:
  const std::string _file;
  const std::string _delimiter;

  // Vector of columns to allow the user to fetch by column index.
  std::vector<std::string> _column_headers;
  // Key: column label.
  // Value: entire column.
  std::unordered_map<std::string, std::vector<std::string>> _csv_data;
}; // class SimpleCSVReader
