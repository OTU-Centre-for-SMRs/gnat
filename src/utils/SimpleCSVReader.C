#include "SimpleCSVReader.h"

#include <fstream>

SimpleCSVReader::SimpleCSVReader(const std::string & file, const std::string & delimiter)
  : _file(file), _delimiter(delimiter)
{
}

bool
SimpleCSVReader::read()
{
  std::ifstream file(_file);
  std::string line;
  if (file.is_open())
  {
    std::size_t current_delim_pos = std::string::npos;
    std::size_t previous_delim_pos = 0u;
    // First line of a CSV file is the header. Parse the header and emplace new vectors in the data
    // map.
    std::getline(file, line);
    current_delim_pos = line.find(_delimiter);
    if (current_delim_pos == std::string::npos && line.size() > 0u)
    {
      // Handle single column edge case.
      _column_headers.emplace_back(line);
      _csv_data.emplace(line, std::vector<std::string>());

      // Grab the rows that contain data.
      while (std::getline(file, line))
        _csv_data[_column_headers[0u]].emplace_back(line);
    }
    else if (current_delim_pos != std::string::npos && line.size() > 0u)
    {
      // Multiple column case.
      std::size_t offset = 0u;
      while (current_delim_pos != std::string::npos)
      {
        _column_headers.emplace_back(line.substr(
            previous_delim_pos + offset, current_delim_pos - (previous_delim_pos + offset)));
        _csv_data.emplace(line.substr(previous_delim_pos + offset,
                                      current_delim_pos - (previous_delim_pos + offset)),
                          std::vector<std::string>());

        previous_delim_pos = current_delim_pos;
        current_delim_pos = line.find(_delimiter, previous_delim_pos + 1);
        offset = 1u;
      }

      _column_headers.emplace_back(line.substr(previous_delim_pos + 1u));
      _csv_data.emplace(line.substr(previous_delim_pos + 1u), std::vector<std::string>());

      // Grab the rows that contain data.
      std::size_t column = 0u;
      while (std::getline(file, line))
      {
        current_delim_pos = line.find(_delimiter);
        previous_delim_pos = 0u;
        column = 0u;
        offset = 0u;
        while (current_delim_pos != std::string::npos)
        {
          _csv_data[_column_headers[column]].emplace_back(line.substr(
              previous_delim_pos + offset, current_delim_pos - (previous_delim_pos + offset)));

          column++;
          previous_delim_pos = current_delim_pos;
          current_delim_pos = line.find(_delimiter, current_delim_pos + 1);
          offset = 1u;
        }

        _csv_data[_column_headers[column]].emplace_back(line.substr(previous_delim_pos));
      }
    }
    else
    {
      // Not a CSV file.
      file.close();
      return false;
    }

    file.close();
    return true;
  }
  else
    return false;
}
