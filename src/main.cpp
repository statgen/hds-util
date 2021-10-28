/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "merge_writer.hpp"

#include <getopt.h>
#include <cstdlib>
#include <cmath>

std::vector<std::string> split_string_to_vector(const char* in, char delim)
{
  std::vector<std::string> ret;
  const char* d = nullptr;
  std::string token;
  const char* s = in;
  const char*const e = in + strlen(in);
  while ((d = std::find(s, e,  delim)) != e)
  {
    ret.emplace_back(std::string(s, d));
    s = d ? d + 1 : d;
  }
  ret.emplace_back(std::string(s,d));
  return ret;
}

class prog_args
{
private:
  std::vector<option> long_options_;
  std::vector<std::string> format_fields_;
  std::vector<std::string> input_paths_;
  std::string output_path_ = "/dev/stdout";
  savvy::file::format output_format_ = savvy::file::format::sav;
  int compression_level_ = 6;
  float min_r2_ = -1.f;
  bool help_ = false;
  bool version_ = false;
public:
  prog_args() :
    long_options_(
      {
        {"format", required_argument, 0, 'f'},
        {"help", no_argument, 0, 'h'},
        {"min-r2", required_argument, 0, 'm'},
        {"output", required_argument, 0, 'o'},
        {"output-format", required_argument, 0, 'O'},
        {"version", no_argument, 0, 'v'},
        {0, 0, 0, 0}
      })
  {
  }

  const std::vector<std::string>& format_fields() const { return format_fields_; }
  const std::vector<std::string>& input_paths() const { return input_paths_; }
  const std::string& output_path() const { return output_path_; }
  savvy::file::format output_format() const { return output_format_; }
  int compression_level() const { return compression_level_; }
  float min_r2() const { return min_r2_; }
  bool help_is_set() const { return help_; }
  bool version_is_set() const { return version_; }

  void print_usage(std::ostream& os)
  {
    os << "Usage: hds-util [opts ...] [input_files.sav ...] \n";
    os << "\n";
    os << " -f, --format         Comma-separated list of FORMAT fields to export (GT, HDS, DS, GP, or SD)\n";
    os << " -h, --help           Print usage\n";
    os << " -m, --min-r2         Minimum r-square threshold\n";
    os << " -o, --output         Output path (default: /dev/stdout)\n";
    os << " -O, --output-format  Output file format (vcf, vcf.gz, bcf, ubcf, sav, usav; default: vcf)\n";
    os << " -v, --version        Print version\n";
    os << std::flush;
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "f:hm:o:O:v", long_options_.data(), &long_index)) != -1)
    {
      char copt = char(opt & 0xFF);
      switch (copt)
      {
      case 'f':
      {
        std::set<std::string> allowed = {"GT", "HDS", "DS", "GP", "SD"};
        format_fields_ = split_string_to_vector(optarg ? optarg : "", ',');
        for (const auto& f : format_fields_)
        {
          if (allowed.find(f) == allowed.end())
            return std::cerr << "Invalid --format value (" << f << ")\n", false;
        }
        break;
      }
      case 'h':
        help_ = true;
        return true;
      case 'm':
        min_r2_ = std::atof(optarg ? optarg : "");
        break;
      case 'o':
        output_path_ = optarg ? optarg : "";
        break;
      case 'O':
      {
        using fmt = savvy::file::format;
        std::string ot = optarg ? optarg : "";
        if (ot == "vcf")
        {
          output_format_ = fmt::vcf;
          compression_level_ = 0;
        }
        else if (ot == "vcf.gz")
        {
          output_format_ = fmt::vcf;
        }
        else if (ot == "bcf")
        {
          output_format_ = fmt::bcf;
        }
        else if (ot == "ubcf")
        {
          output_format_ = fmt::bcf;
          compression_level_ = 0;
        }
        else if (ot == "sav")
        {
          output_format_ = fmt::sav;
        }
        else if (ot == "usav")
        {
          output_format_ = fmt::sav;
          compression_level_ = 0;
        }
        else
        {
          std::cerr << "Invalid --output-format: " << ot << std::endl;
          return false;
        }
        break;
      }
      case 'v':
        version_ = true;
        return true;
      default:
        return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count == 0)
    {
      input_paths_.push_back("/dev/stdin");
    }
    else
    {
      input_paths_.reserve(remaining_arg_count);
      for (std::size_t i = 0; i < remaining_arg_count; ++i)
        input_paths_.push_back(argv[optind + i]);
    }

    return true;
  }
};



bool is_empirical_file(const savvy::reader& rdr)
{
  for (auto it = rdr.format_headers().begin(); it != rdr.format_headers().end(); ++it)
  {
    if (it->id == "LDS")
      return true;
  }
  return false;
}

std::vector<std::string> format_fields_in_headers(const savvy::reader& rdr)
{
  std::vector<std::string> ret(rdr.format_headers().size());
  std::size_t i = 0;
  for (auto it = rdr.format_headers().begin(); it != rdr.format_headers().end(); ++it)
    ret[i++] = it->id;
  return ret;
}

int main(int argc, char** argv)
{
  prog_args args;
  if (!args.parse(argc, argv))
  {
    args.print_usage(std::cerr);
    return EXIT_FAILURE;
  }

  if (args.help_is_set())
  {
    args.print_usage(std::cout);
    return EXIT_SUCCESS;
  }

  if (args.version_is_set())
  {
    std::cout << "hds-util v" << VERSION << std::endl;
    return EXIT_SUCCESS;
  }

  std::list<savvy::reader> input_files;
  std::vector<std::vector<std::string>> sample_ids(args.input_paths().size());
  std::vector<bool> empirical(args.input_paths().size());
  for (std::size_t i = 0; i < args.input_paths().size(); ++i)
  {
    input_files.emplace_back(args.input_paths()[i]);
    sample_ids[i] = input_files.back().samples();
    empirical[i] = is_empirical_file(input_files.back());
  }

  auto format_fields = args.format_fields();

  std::size_t emp_cnt = std::count(empirical.begin(), empirical.end(), true);
  if (emp_cnt && emp_cnt != empirical.size())
    return std::cerr << "Error: cannot combine empirical files with non-empirical files\n", EXIT_FAILURE;

  assert(input_files.size());
  if (emp_cnt)
    format_fields = {"GT", "LDS"};
  else if (format_fields.empty())
    format_fields = format_fields_in_headers(input_files.front());

  if (!format_fields.empty() && emp_cnt)
    return std::cerr << "Error: --format option not valid with empirical files\n", EXIT_FAILURE;

  merge_writer output(args.output_path(),
    args.output_format(), args.compression_level(),
    input_files.front().headers(),
    sample_ids,
    format_fields,
    args.min_r2());

  return output.merge_files(input_files) ? EXIT_SUCCESS : EXIT_FAILURE;
}
