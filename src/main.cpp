/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <savvy/reader.hpp>
#include <savvy/writer.hpp>

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
  static const int default_compression_level = 3;
  static const int default_block_size = 2048;

  std::vector<option> long_options_;
  std::vector<std::string> generate_fields_;
  std::string input_path_;
  std::string output_path_;
  savvy::file::format output_format_;
  int compression_level_ = 0;
  bool help_ = false;
public:
  prog_args() :
    long_options_(
      {
        {"generate", required_argument, 0, 'g'},
        {"help", no_argument, 0, 'h'},
        {"output", required_argument, 0, 'o'},
        {"output-format", required_argument, 0, 'O'},
        {0, 0, 0, 0}
      })
  {
  }

  const std::vector<std::string>& generate_fields() const { return generate_fields_; }
  const std::string& input_path() const { return input_path_; }
  const std::string& output_path() const { return output_path_; }
  savvy::file::format output_format() const { return output_format_; }
  int compression_level() const { return compression_level_; }
  bool help_is_set() const { return help_; }

  void print_usage(std::ostream& os)
  {
    os << "Usage: sav index [opts ...] <in.sav> \n";
    os << "\n";
    os << " -g, --generate       Comma-separated list of FORMAT fields to generate (GT, DS, GP, or SD)\n";
    os << " -h, --help           Print usage\n";
    os << " -o, --output         Output path (default: /dev/stdout)\n";
    os << " -O, --output-format  Output format (default: vcf)\n";
    os << std::flush;
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "g:ho:O:", long_options_.data(), &long_index)) != -1)
    {
      char copt = char(opt & 0xFF);
      switch (copt)
      {
      case 'g':
      {
        std::set<std::string> allowed = {"GT", "DS", "GP", "SD"};
        generate_fields_ = split_string_to_vector(optarg ? optarg : "", ',');
        for (const auto& g : generate_fields_)
        {
          if (allowed.find(g) == allowed.end())
            return std::cerr << "Invalid --generate value (" << g << ")\n", false;
        }
        break;
      }
      case 'h':
        help_ = true;
        return true;
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
        }
        else if (ot == "vcf.gz")
        {
          output_format_ = fmt::vcf;
          compression_level_ = 6;
        }
        else if (ot == "bcf")
        {
          output_format_ = fmt::bcf;
          compression_level_ = 6;
        }
        else if (ot == "ubcf")
        {
          output_format_ = fmt::bcf;
        }
        else if (ot == "ubcf")
        {
          output_format_ = fmt::sav;
          compression_level_ = 6;
        }
        else
        {
          std::cerr << "Invalid --output-format: " << ot << std::endl;
          return false;
        }
        break;
      }
      default:
        return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count == 1)
    {
      input_path_ = argv[optind];
    }
    else if (remaining_arg_count < 1)
    {
      std::cerr << "Too few arguments\n";
      return false;
    }
    else
    {
      std::cerr << "Too many arguments\n";
      return false;
    }

    if (input_path_ == "/dev/stdin" || input_path_ == "/dev/fd/0")
    {
      std::cerr << "Input SAV file cannot be stdin\n";
      return false;
    }

    return true;
  }
};

class field_generator
{
public:
  field_generator(const std::vector<std::string>& fields, std::size_t n_samples) :
    fields_(fields),
    n_samples_(n_samples)
  {
  }

  bool operator()(savvy::variant& record)
  {
    if (!record.get_format("HDS", hds_vec_))
      return false;

    std::size_t stride = hds_vec_.size() / n_samples_;

    for (const auto& field : fields_)
    {
      if (field == "GT")
      {
        gt_vec_.clear();
        gt_vec_.resize(hds_vec_.size());
        for (std::size_t i = 0; i < hds_vec_.size(); ++i)
        {
          if (hds_vec_[i] == 0.f)
            continue;
          gt_vec_[i] = std::isnan(hds_vec_[i]) ? hds_vec_[i] : std::int8_t(std::round(hds_vec_[i]));
        }

        record.set_format("GT", gt_vec_);
      }
      else if (field == "DS")
      {
        ds_vec_.clear();
        ds_vec_.resize(n_samples_);
        for (std::size_t i = 0; i < n_samples_; ++i)
        {
          float ds = hds_vec_[i * stride];
          for (std::size_t j = 1; j < stride; ++j)
          {
            float v = hds_vec_[i * stride + j];
            if (savvy::typed_value::is_end_of_vector(v))
              break;
            ds += v;
          }

          if (std::isnan(ds) || ds)
            ds_vec_[i] = ds;
        }
        record.set_format("DS", ds_vec_);
      }
      else if (field == "GP")
      {
        if (stride == 1)
        {
          // All samples are haploid
          dense_float_vec_.resize(n_samples_ * 2);
          for (std::size_t i = 0; i < n_samples_; ++i)
          {
            std::size_t dest_idx = i * 2;
            dense_float_vec_[dest_idx] = 1.f - hds_vec_[i];
            dense_float_vec_[dest_idx + 1] = hds_vec_[i];
          }
        }
        else if (stride == 2)
        {
          dense_float_vec_.resize(n_samples_ * 3);
          for (std::size_t i = 0; i < n_samples_; ++i)
          {
            std::size_t src_idx = i * 2;
            std::size_t dest_idx = i * 3;
            float x = hds_vec_[src_idx];
            float y = hds_vec_[src_idx + 1];
            if (savvy::typed_value::is_end_of_vector(y))
            {
              // haploid
              dense_float_vec_[dest_idx] = 1.f - x;
              dense_float_vec_[dest_idx + 1] = x;
              dense_float_vec_[dest_idx + 2] = y;
            }
            else
            {
              // diploid
              dense_float_vec_[dest_idx] = (1.f - x) * (1.f - y);
              dense_float_vec_[dest_idx + 1] = x * (1.f - y) + y * (1.f - x);
              dense_float_vec_[dest_idx + 2] =  x * y;
            }
          }
        }
        else
        {
          std::cerr << "Error: only haploid and diploid samples are supported when generating GP\n";
          return false;
        }
        record.set_format("GP", dense_float_vec_);
      }
      else if (field == "SD")
      {
        dense_float_vec_.resize(n_samples_);
        if (stride == 1)
        {
          // All samples are haploid
          for (std::size_t i = 0; i < n_samples_; ++i)
          {
            dense_float_vec_[i] = hds_vec_[i] * (1.f - hds_vec_[i]);
          }
        }
        else if (stride == 2)
        {
          for (std::size_t i = 0; i < hds_vec_.size(); i+=2)
          {
            float x = hds_vec_[i];
            float y = hds_vec_[i + 1];
            if (savvy::typed_value::is_end_of_vector(y)) // haploid
              dense_float_vec_[i / 2] = x * (1.f - x);
            else // diploid
              dense_float_vec_[i / 2] = x * (1.f - x) + y * (1.f - y);
          }
        }
        else
        {
          std::cerr << "Error: only haploid and diploid samples are supported when generating SD\n";
          return false;
        }
        record.set_format("SD", dense_float_vec_);
      }
      else
      {
        std::cerr << "Error: " << field << " field not supported\n";
        return false;
      }
    }
    return true;
  }
private:
  std::size_t n_samples_ = 0;
  std::vector<std::string> fields_;
  std::vector<float> hds_vec_;
  std::vector<float> dense_float_vec_;
  savvy::compressed_vector<std::int8_t> gt_vec_;
  savvy::compressed_vector<float> ds_vec_;
};

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

  std::unordered_map<std::string, std::string> header_map = {
    {"GT", "<ID=GT,Number=1,Type=String,Description=\"Genotype\">"},
    {"DS", "<ID=DS,Number=1,Type=Float,Description=\"Estimated Alternate Allele Dosage : [P(0/1)+2*P(1/1)]\">"},
    {"GP", "<ID=GP,Number=G,Type=Float,Description=\"Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1 \">"},
    {"SD", "<ID=SD,Number=1,Type=Float,Description=\"Variance of Posterior Genotype Probabilities\">"}
  };

  savvy::reader input(args.input_path());

  auto hdrs = input.headers();

  std::set<std::string> existing_format_headers;
  for (const auto& h : input.format_headers())
    existing_format_headers.insert(h.id);

  for (auto g : args.generate_fields())
  {
    if (existing_format_headers.find(g) == existing_format_headers.end())
      hdrs.emplace_back("FORMAT", header_map[g]);
  }

  savvy::writer output(args.output_path(), args.output_format(), hdrs, input.samples(), args.compression_level());

  savvy::variant record;
  field_generator generate_fields(args.generate_fields(), input.samples().size());
  savvy::typed_value tv;

  while (input >> record && output)
  {
    generate_fields(record);

    if (args.output_format() == savvy::file::format::sav)
    {
      // ensure HDS is sparse
      for (auto it = record.format_fields().begin(); it != record.format_fields().end(); ++it)
      {
        if (it->first == "HDS")
        {
          if (!it->second.is_sparse())
          {
            it->second.copy_as_sparse(tv);
            record.set_format("HDS", tv);
          }
          break;
        }
      }
    }

    output << record;
  }

  if (!output)
    return std::cerr << "Error: write failure\n", EXIT_FAILURE;

  if (input.bad())
    return std::cerr << "Error: read failure\n", EXIT_FAILURE;

  return EXIT_SUCCESS;
}
