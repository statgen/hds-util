/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef HDSUTIL_MERGER_WRITER_HPP
#define HDSUTIL_MERGER_WRITER_HPP

#include <savvy/writer.hpp>
#include <savvy/reader.hpp>

class merge_writer
{
private:
  savvy::writer out_file_;
  std::unordered_set<std::string> fmt_field_set_;
  std::size_t n_samples_ = 0;
  float min_r2_ = -1.f;

  // buffers
  savvy::compressed_vector<std::int8_t> sparse_gt_;
  std::vector<float> dense_float_vec_;
  std::vector<float> dense_zero_vec_;
public:
  typedef std::vector<std::string> str_vec;
  merge_writer(const std::string& file_path,
    savvy::file::format file_format,
    std::uint8_t out_compression,
    const std::vector<std::pair<std::string, std::string>>& reader_headers,
    const std::vector<std::vector<std::string>>& sample_ids,
    const std::vector<std::string>& fmt_fields,
    float min_r2 = -1.f)
    :
    out_file_(file_path, file_format, update_headers(fmt_fields, reader_headers, sample_ids.size()),
      std::accumulate(sample_ids.begin(), sample_ids.end(), str_vec(), [](str_vec& l, const str_vec& r) { l.insert(l.end(), r.begin(), r.end()); return l; }),
      out_compression),
    fmt_field_set_(fmt_fields.begin(), fmt_fields.end()),
    n_samples_(std::accumulate(sample_ids.begin(), sample_ids.end(), std::size_t(0), [](std::size_t& l, const str_vec& r) { return l + r.size(); })),
    min_r2_(min_r2)
  {

  }

  static std::vector<std::pair<std::string, std::string>> update_headers(const std::vector<std::string>& fmt_fields, const std::vector<std::pair<std::string, std::string>>& reader_headers, std::size_t n_source_files)
  {
    std::unordered_set<std::string> known_info_fields = {"AF","MAF","R2","ER2","AVG_CS"};
    std::vector<std::pair<std::string, std::string>> known_info_headers = {
      {"INFO", "<ID=AF,Number=1,Type=Float,Description=\"Estimated Alternate Allele Frequency\">"},
      {"INFO", "<ID=MAF,Number=1,Type=Float,Description=\"Estimated Minor Allele Frequency\">"},
      {"INFO", "<ID=AVG_CS,Number=1,Type=Float,Description=\"Average Call Score\">"},
      {"INFO", "<ID=R2,Number=1,Type=Float,Description=\"Estimated Imputation Accuracy (R-square)\">"},
      {"INFO", "<ID=ER2,Number=" + std::to_string(n_source_files) + ",Type=Float,Description=\"Empirical (Leave-One-Out) R-square (one value per imputation group)\">"}};

    std::time_t t = std::time(nullptr);
    char datestr[11];
    std::string filedate(datestr, std::strftime(datestr, sizeof(datestr), "%Y%m%d", std::localtime(&t)));
    assert(filedate.size());

    std::vector<std::pair<std::string, std::string>> ret;
    ret.reserve(reader_headers.size() + fmt_fields.size());

    for (auto it = reader_headers.begin(); it != reader_headers.end(); ++it)
    {
      if (it->first == "filedate")
      {
        ret.emplace_back(it->first, filedate);
      }
      else if (it->first == "FORMAT")
      {
        // FORMAT headers are added later
      }
      else if (it->first == "phasing")
      {
        // overwrite phasing header later
      }
      else if (it->first == "INFO")
      {
        if (known_info_fields.find(savvy::parse_header_sub_field(it->second, "ID")) == known_info_fields.end())
          ret.push_back(*it);
      }
      else
      {
        ret.push_back(*it);
      }
    }

    ret.insert(ret.end(), known_info_headers.begin(), known_info_headers.end());

    for (const auto& f : fmt_fields)
    {
      if (f == "GT")
        ret.emplace_back("FORMAT", "<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
      else if (f == "DS")
        ret.emplace_back("FORMAT", "<ID=DS,Number=1,Type=Float,Description=\"Estimated Alternate Allele Dosage : [P(0/1)+2*P(1/1)]\">");
      else if (f == "HDS")
        ret.emplace_back("FORMAT", "<ID=HDS,Number=.,Type=Float,Description=\"Estimated Haploid Alternate Allele Dosage \">");
      else if (f == "GP")
        ret.emplace_back("FORMAT", "<ID=GP,Number=G,Type=Float,Description=\"Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1 \">");
      else if (f == "SD")
        ret.emplace_back("FORMAT", "<ID=SD,Number=1,Type=Float,Description=\"Variance of Posterior Genotype Probabilities\">");
      else if (f == "LDS")
        ret.emplace_back("FORMAT", "<ID=LDS,Number=.,Type=Float,Description=\"Leave-one-out Imputed Dosage : Estimated Haploid Alternate Allele Dosage assuming site was NOT genotyped\">");
    }

    ret.emplace_back("phasing", "full");

    // TODO: command string

    return ret;
  }

  bool merge_files(std::list<savvy::reader>& temp_files)
  {
    if (temp_files.empty())
      return false;

    savvy::variant out_var;
    savvy::compressed_vector<float> pasted_hds;
    std::vector<savvy::compressed_vector<float>> partial_hds_vecs(temp_files.size());
    std::vector<float> partial_hds_dense;

    std::vector<int8_t> pasted_gt;
    std::vector<savvy::compressed_vector<std::int8_t>> partial_gt_vecs(temp_files.size());

    std::vector<float> pasted_er2;

    bool is_empirical = fmt_field_set_.find("LDS") != fmt_field_set_.end();

    int good_count = temp_files.size();
    while (good_count == temp_files.size() && out_file_)
    {
      pasted_hds.clear();
      pasted_gt.clear();
      pasted_er2.clear();

      std::size_t max_ploidy = 0;
      std::vector<std::size_t> ploidies(temp_files.size());
      std::uint32_t prev_pos;
      std::vector<std::string> prev_alts;
      good_count = 0;
      std::size_t i = 0;
      for (auto it = temp_files.begin(); it != temp_files.end(); ++it,++i)
      {
        if (it->read(out_var))
        {
          if (i == 0)
          {
            prev_pos = out_var.pos();
            prev_alts = out_var.alts();
          }
          else if (prev_pos != out_var.pos() || prev_alts != out_var.alts())
          {
            return std::cerr << "Error: site lists do not match\n", false;
          }

          ++good_count;

          float er2;
          if (out_var.get_info("ER2", er2))
            pasted_er2.push_back(er2);

          if ((!is_empirical && !out_var.get_format("HDS", partial_hds_vecs[i])) || (is_empirical && !out_var.get_format("LDS", partial_hds_vecs[i])))
            return std::cerr << "Error: " << (is_empirical ? "LDS" : "HDS") << " must be present\n", false;
          ploidies[i] = partial_hds_vecs[i].size() / it->samples().size();
          max_ploidy = std::max(max_ploidy, ploidies[i]);

          if (is_empirical)
          {
            if (!out_var.get_format("GT", partial_gt_vecs[i]))
              return std::cerr << "Error: GT must be present in empirical files\n", false;

            if (partial_hds_vecs[i].size() != partial_gt_vecs[i].size())
              return std::cerr << "Error: Length of GT does not match length of LDS\n", false;

            if (ploidies[i] != partial_gt_vecs[i].size() / it->samples().size())
              return std::cerr << "Error: Ploidy of GT does not match ploidy of LDS\n", false;
          }
        }
      }

      if (good_count)
      {
        if (good_count < temp_files.size())
          return std::cerr << "Error: record mismatch in temp files\n", false;

        pasted_hds.resize(max_ploidy * n_samples_);
        if (is_empirical)
          pasted_gt.resize(max_ploidy * n_samples_);
        std::size_t off = 0;
        for (std::size_t i = 0; i < partial_hds_vecs.size(); ++i)
        {
          if (ploidies[i] < max_ploidy)
          {
            partial_hds_dense.clear();
            partial_hds_dense.resize(partial_hds_vecs[i].size() / ploidies[i] * max_ploidy);
            for (std::size_t j = 0; j < partial_hds_vecs[i].size() / ploidies[i]; ++j)
            {
              for (std::size_t k = ploidies[i]; k < max_ploidy; ++k)
                partial_hds_dense[j * max_ploidy + k] = savvy::typed_value::end_of_vector_value<float>();
            }

            for (auto jt = partial_hds_vecs[i].begin(); jt != partial_hds_vecs[i].end(); ++jt)
              partial_hds_dense[jt.offset() / ploidies[i] * max_ploidy + jt.offset() % ploidies[i]] = *jt;

            for (std::size_t j = 0; j < partial_hds_dense.size(); ++j)
            {
              if (partial_hds_dense[j] == 0.f) continue;
              pasted_hds[off + j] = partial_hds_dense[j];
            }

            if (is_empirical)
            {
              for (std::size_t j = 0; j < partial_gt_vecs[i].size() / ploidies[i]; ++j)
              {
                for (std::size_t k = ploidies[i]; k < max_ploidy; ++k)
                  pasted_gt[off + j * max_ploidy + k] = savvy::typed_value::end_of_vector_value<std::int8_t>();
              }

              for (auto jt = partial_gt_vecs[i].begin(); jt != partial_gt_vecs[i].end(); ++jt)
                pasted_gt[off + jt.offset() / ploidies[i] * max_ploidy + jt.offset() % ploidies[i]] = *jt;
            }

            off += partial_hds_dense.size();
          }
          else
          {
            for (auto jt = partial_hds_vecs[i].begin(); jt != partial_hds_vecs[i].end(); ++jt)
              pasted_hds[off + jt.offset()] = *jt;

            if (is_empirical)
            {
              for (auto jt = partial_gt_vecs[i].begin(); jt != partial_gt_vecs[i].end(); ++jt)
                pasted_gt[off + jt.offset()] = *jt;
            }

            off += partial_hds_vecs[i].size();
          }
        }

        if (!is_empirical)
          set_info_fields(out_var, pasted_hds);

        out_var.remove_info("S_CS");
        out_var.remove_info("S_X");
        out_var.remove_info("S_XX");
        out_var.remove_info("LOO_S_X");
        out_var.remove_info("LOO_S_XX");
        out_var.remove_info("LOO_S_Y");
        out_var.remove_info("LOO_S_YY");
        out_var.remove_info("LOO_S_XY");
        out_var.remove_info("AN");

        if (has_good_r2(out_var))
        {
          if (pasted_er2.size())
            out_var.set_info("ER2", pasted_er2);

          if (is_empirical)
          {
            out_var.set_format("GT", pasted_gt);
            out_var.set_format("LDS", pasted_hds);
          }
          else
          {
            set_format_fields(out_var, pasted_hds);
          }

          // Remove unexpected (unmerged) FORMAT fields
          std::vector<std::string> fmt_fields_to_remove;
          for (auto fmt = out_var.format_fields().begin(); fmt != out_var.format_fields().end(); ++fmt)
          {
            if (fmt_field_set_.find(fmt->first) == fmt_field_set_.end())
              fmt_fields_to_remove.push_back(fmt->first);
          }

          for (auto fmt = fmt_fields_to_remove.begin(); fmt != fmt_fields_to_remove.end(); ++fmt)
            out_var.set_format(*fmt, {});

          out_file_ << out_var;
        }
      }
    }

    int bad_count = 0;
    for (auto it = temp_files.begin(); it != temp_files.end(); ++it)
      bad_count += (int)it->bad();

    if (bad_count || !out_file_.good())
      return std::cerr << "Error: I/O failed while merging" << std::endl, false;

    return true;
  }
private:
  bool has_good_r2(savvy::site_info& site)
  {
    if (min_r2_ >= 0.f)
    {
      float r2 = -1.f;
      site.get_info("R2", r2);
      if (r2 >= min_r2_)
        return true;
      return false;
    }

    return true;
  }

  void set_r2_info_field(savvy::variant& out_var, float s_x, float s_xx, std::size_t n)
  {
    float af = s_x / float(n);
    float r2 = savvy::typed_value::missing_value<float>();
    if (af > 0.f && af < 1.f)
      r2 = ((s_xx - s_x * s_x / n) / n) / (af * (1.f - af));

    out_var.set_info("R2", r2);
  }

//  void set_er2_info_field(savvy::variant& out_var, float s_x, float s_xx, float s_y, float s_yy, float s_xy, std::size_t n)
//  {
//    //                         n * Sum xy - Sum x * Sum y
//    //  r = -------------------------------------------------------------------
//    //      Sqrt(n * Sum xx - Sum x * Sum x) * Sqrt(n * Sum yy - Sum y * Sum y)
//    float emp_r = (n * s_xy - s_x * s_y) / (std::sqrt(n * s_xx - s_x * s_x) * std::sqrt(n * s_yy - s_y * s_y));
//    out_var.set_info("ER2", std::isnan(emp_r) ? savvy::typed_value::missing_value<float>() : emp_r * emp_r);
//  }

  struct plus_ignore_missing
  {
    float operator()(const float& l, const float& r)
    {
      if (std::isnan(r))
        return l;
      if (std::isnan(l))
        return r;
      return l + r;
    }
  };

  void set_info_fields(savvy::variant& out_var, const savvy::compressed_vector<float>& sparse_dosages)
  {
    std::size_t n = sparse_dosages.size();
    assert(n);

    float s_x = std::accumulate(sparse_dosages.begin(), sparse_dosages.end(), 0.f, plus_ignore_missing());
    float s_xx = std::inner_product(sparse_dosages.begin(), sparse_dosages.end(), sparse_dosages.begin(), 0.f, plus_ignore_missing(), std::multiplies<float>());
    float s_cs(sparse_dosages.size() - sparse_dosages.non_zero_size());
    for (auto it = sparse_dosages.begin(); it != sparse_dosages.end(); ++it)
    {
      if (savvy::typed_value::is_end_of_vector(*it))
        --n;
      else
        s_cs += *it > 0.5f ? *it : 1.f - *it;
    }

    float af = s_x / n;

    out_var.set_info("AF", af);
    out_var.set_info("MAF", af > 0.5f ? 1.f - af : af);
    out_var.set_info("AVG_CS", s_cs / n);
    set_r2_info_field(out_var, s_x, s_xx, n);
  }

  void set_format_fields(savvy::variant& out_var, savvy::compressed_vector<float>& sparse_dosages)
  {
    std::size_t stride = sparse_dosages.size() / n_samples_;

    if (fmt_field_set_.find("GT") != fmt_field_set_.end())
    {
      if (out_var.format_fields().size() && out_var.format_fields().front().first != "GT")
      {
        // Ensure GT is first FORMAT field.
        std::vector<std::string> fmt_keys;
        fmt_keys.reserve(out_var.format_fields().size());
        for (auto it = out_var.format_fields().begin(); it != out_var.format_fields().end(); ++it)
          fmt_keys.push_back(it->first);

        for (auto it = fmt_keys.begin(); it != fmt_keys.end(); ++it)
          out_var.set_format(*it, {});
      }

      sparse_gt_.assign(sparse_dosages.value_data(), sparse_dosages.value_data() + sparse_dosages.non_zero_size(), sparse_dosages.index_data(), sparse_dosages.size(), [](float v)
      {
        if (savvy::typed_value::is_end_of_vector(v))
          return savvy::typed_value::end_of_vector_value<std::int8_t>();
        return std::int8_t(v < 0.5f ? 0 : 1);
      });
      out_var.set_format("GT", sparse_gt_);
    }

    if (fmt_field_set_.find("HDS") != fmt_field_set_.end())
    {
      out_var.set_format("HDS", sparse_dosages);
    }
    else
    {
      out_var.set_format("HDS", {});
    }

    if (fmt_field_set_.find("GP") != fmt_field_set_.end() || fmt_field_set_.find("SD") != fmt_field_set_.end())
    {
      // set dense dosage vector
      dense_zero_vec_.resize(sparse_dosages.size());
      for (auto it = sparse_dosages.begin(); it != sparse_dosages.end(); ++it)
        dense_zero_vec_[it.offset()] = *it;

      std::vector<float>& dense_hds = dense_zero_vec_;

      if (fmt_field_set_.find("GP") != fmt_field_set_.end())
      {
        if (stride == 1)
        {
          // All samples are haploid
          dense_float_vec_.resize(n_samples_ * 2);
          for (std::size_t i = 0; i < n_samples_; ++i)
          {
            std::size_t dest_idx = i * 2;
            dense_float_vec_[dest_idx] = 1.f - dense_hds[i];
            dense_float_vec_[dest_idx + 1] = dense_hds[i];
          }
        }
        else if (stride == 2)
        {
          dense_float_vec_.resize(n_samples_ * 3);
          for (std::size_t i = 0; i < n_samples_; ++i)
          {
            std::size_t src_idx = i * 2;
            std::size_t dest_idx = i * 3;
            float x = dense_hds[src_idx];
            float y = dense_hds[src_idx + 1];
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
              dense_float_vec_[dest_idx + 2] = x * y;
            }
          }
        }

        out_var.set_format("GP", dense_float_vec_);
      }
      else
      {
        out_var.set_format("GP", {});
      }

      if (fmt_field_set_.find("SD") != fmt_field_set_.end())
      {
        dense_float_vec_.resize(n_samples_);
        if (stride == 1)
        {
          // All samples are haploid
          for (std::size_t i = 0; i < n_samples_; ++i)
          {
            dense_float_vec_[i] = dense_hds[i] * (1.f - dense_hds[i]);
          }

          out_var.set_format("SD", dense_float_vec_);
        }
        else if (stride == 2)
        {
          for (std::size_t i = 0; i < dense_hds.size(); i += 2)
          {
            float x = dense_hds[i];
            float y = dense_hds[i + 1];
            if (savvy::typed_value::is_end_of_vector(y)) // haploid
              dense_float_vec_[i / 2] = x * (1.f - x);
            else // diploid
              dense_float_vec_[i / 2] = x * (1.f - x) + y * (1.f - y);
          }

          out_var.set_format("SD", dense_float_vec_);
        }
        else
        {
          // TODO: suppress error excessive error messages
          std::cerr << "Error: only haploid and diploid samples are supported when generating SD\n";
        }
      }
      else
      {
        out_var.set_format("SD", {});
      }

      // unset dense dosage vector
      for (auto it = sparse_dosages.begin(); it != sparse_dosages.end(); ++it)
        dense_zero_vec_[it.offset()] = 0.f;
    }
    else
    {
      out_var.set_format("GP", {});
      out_var.set_format("SD", {});
    }


    if (fmt_field_set_.find("DS") != fmt_field_set_.end())
    {
      savvy::stride_reduce(sparse_dosages, sparse_dosages.size() / n_samples_, savvy::plus_eov<float>());
      out_var.set_format("DS", sparse_dosages);
    }
    else
    {
      out_var.set_format("DS", {});
    }
  }
};

#endif // HDSUTIL_MERGER_WRITER_HPP