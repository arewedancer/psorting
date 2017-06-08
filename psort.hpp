#include <iostream>
#include <algorithm>
#include <iomanip>
#include <unistd.h>
#include <string.h>
#include <cassert>
#include <vector>
#include <memory>
#include <cmath>
#include <cstdint>
// TBB header
#include "tbb/tbb.h"
#include "tbb/task_scheduler_init.h"
#include "merge_sort_tbb.h"
#include "merge_tbb.h"
#include <thread>
#include "sample_sort_util.h"
#include "quicksort_util.h"
#include "quicksort_tbb_task.h"
#include "map_keys_to_bins_tbb.h"
#include "bin_tbb.h"
#include "repack_and_subsort_tbb.h"
#include "sample_sort.h"
#include "tbb/parallel_sort.h"
#include "RankOptions.hpp"

inline bool is_null(double num)
{
	  return (std::signbit(num) && num == 0.0);
};

struct Functor {
	//const std::vector<long>& dim_cols;
	const std::vector<long>& dim_cols;
	const std::vector<long>& measure_cols;
	inline explicit Functor(const std::vector<long>& dim_cols_, const std::vector<long>& measure_cols_):
		dim_cols(dim_cols_), measure_cols(measure_cols_){};

	inline bool operator() (const double* a, const double* b) const { 
		for (int i = 0; i < dim_cols.size(); ++i )
		{
			if (is_null(a[dim_cols[i]]) ^ is_null(b[dim_cols[i]]) )
				return is_null(a[dim_cols[i]]) ? true : false;

			if ( a[dim_cols[i]] == b[dim_cols[i]]) continue;
			//always sort partitions ascending
			return (a[dim_cols[i]] > b[dim_cols[i]]); 
		}

		for (int i = 0; i < measure_cols.size(); ++i)
		{
			if (is_null(a[std::abs(measure_cols[i])]) ^ is_null(b[std::abs(measure_cols[i])]) )
			{
				return is_null(a[std::abs(measure_cols[i])]) ^ (!std::signbit(measure_cols[i]));
			};
			if ( a[std::abs(measure_cols[i])] == b[std::abs(measure_cols[i])] ) continue; 
			return  (a[std::abs(measure_cols[i])]> b[std::abs(measure_cols[i])]) ^ (!std::signbit(measure_cols[i]));
		}
	};
};

struct Pair
{
	long start;long end;
	Pair(long start_, long end_):start(start_),end(end_){};
	Pair(const Pair& p):start(p.start), end(p.end){}; 
	Pair& operator = (const Pair &t) 
	{
		start = t.start; end=t.end;
	};  
};

class Util
{
	private:
		inline static double** get_data(char * data)
		{
			char * d = data + 3 * sizeof(uint32_t);
			return (double**) d;
		};

		inline static void set_row(char * data, uint32_t row)
		{
			uint32_t *d = (uint32_t*)data;
			d[0] = (uint32_t) row;
		};
		inline static bool is_same(const double* a, const double* b, const std::vector<long>& cols)
		{
			for(auto idx: cols)
			{
				if (a[idx] != b[idx])  return false;
			}
			return true;
		};
		typedef uint32_t (*OP_FUNC)(char * result, RankOptions& rankOptions, std::vector<Pair>& index,
				std::vector<long>& measure_values, std::vector<double>& sum);
		static uint32_t mark_rank(char *result, RankOptions& rankOptions, std::vector<Pair>& index, 
				std::vector<long>& measure_values, std::vector<double>& sum);
		static uint32_t mark_denserank(char *result, RankOptions& rankOptions, std::vector<Pair>& index, 
				std::vector<long>& measure_values, std::vector<double>& sum);
		static uint32_t mark_percrank(char *result, RankOptions& rankOptions, std::vector<Pair>& index, 
				std::vector<long>& measure_values, std::vector<double>& sum);
		static uint32_t mark_tb_count(char *result, RankOptions& rankOptions, std::vector<Pair>& index, 
				std::vector<long>& measure_values, std::vector<double>& sum);
		static uint32_t mark_tb_percent(char *result, RankOptions& rankOptions, std::vector<Pair>& index, 
				std::vector<long>& measure_values, std::vector<double>& sum);
		static uint32_t mark_tb_percent_on(char *result, RankOptions& rankOptions, std::vector<Pair>& index, 
				std::vector<long>& measure_values, std::vector<double>& sum);
		static uint32_t mark_ntile(char *result, RankOptions& rankOptions, std::vector<Pair>& index, 
				std::vector<long>& measure_values, std::vector<double>& sum);
		static uint32_t mark_rank_p(char *result, RankOptions& rankOptions, std::vector<Pair>& index, 
				std::vector<long>& measure_values, std::vector<double>& sum);
		static uint32_t mark_denserank_p(char *result, RankOptions& rankOptions, std::vector<Pair>& index, 
				std::vector<long>& measure_values, std::vector<double>& sum);
		static uint32_t mark_percrank_p(char *result, RankOptions& rankOptions, std::vector<Pair>& index, 
				std::vector<long>& measure_values, std::vector<double>& sum);
		static uint32_t mark_ntile_p(char *result, RankOptions& rankOptions, std::vector<Pair>& index, 
				std::vector<long>& measure_values, std::vector<double>& sum);
		static void psum(RankOptions& rankOptions, std::vector<Pair>& index);
		static void psum(RankOptions& rankOptions, std::vector<double>& sum, std::vector<Pair>& index);

		static constexpr OP_FUNC op_map_parallel[7]  = {mark_tb_count, mark_tb_percent, mark_tb_percent_on, mark_rank_p, mark_denserank_p, mark_percrank_p, mark_ntile_p};
		static constexpr OP_FUNC op_map[7]  = {mark_tb_count, mark_tb_percent, mark_tb_percent_on, mark_rank, mark_denserank, mark_percrank, mark_ntile};
	public:
		static char* new_new(const long row, const long measure_col_cnt, const long dim_col_cnt);
		static char* sort(RankOptions& rankOptions);
};

inline void copy_results(RankOptions rankOptions, double current_type, long idx, double ** result)
{
	for (long j = 0; j < rankOptions.result_cols.size(); ++j)
	{
		if ( rankOptions.result_cols[j] != Rank_Result_Index )
			result[idx][j] = rankOptions.get_data()[idx][rankOptions.result_cols[j]];
		else
			result[idx][j] = current_type;
	};
};

uint32_t Util::mark_rank(char *result, RankOptions& rankOptions, std::vector<Pair>& index, 
		std::vector<long>& measure_values, std::vector<double>& sum)
{
	for(long i = 0; i < index.size(); ++i)
	{
		long pre_idx = index[i].start; long current_type = 1;
		long idx = index[i].start;
		while (idx <= index[i].end)
		{
			if (!is_same(rankOptions.get_data()[pre_idx], rankOptions.get_data()[idx], measure_values))
			{
				pre_idx = idx;
				current_type = idx - index[i].start + 1;
			}
			copy_results(rankOptions, current_type, idx, get_data(result));
			idx++;
		};
	}
	return rankOptions.get_row();
};

uint32_t Util::mark_rank_p(char *result, RankOptions& rankOptions, std::vector<Pair>& index, 
		std::vector<long>& measure_values, std::vector<double>& sum)
{
	tbb::parallel_for(size_t(0), index.size(), size_t(1) , [&](size_t i) {
			long pre_idx = index[i].start; long current_type = 1;
			long idx = index[i].start;
			while (idx <= index[i].end)
			{
			if (!is_same(rankOptions.get_data()[pre_idx], rankOptions.get_data()[idx], measure_values))
			{
			pre_idx = idx;
			current_type = idx - index[i].start + 1;
			}
			copy_results(rankOptions, current_type, idx, get_data(result));
			idx++;
			};
			});
	return rankOptions.get_row();
};
uint32_t Util::mark_percrank(char *result, RankOptions& rankOptions, std::vector<Pair>& index, 
		std::vector<long>& measure_values, std::vector<double>& sum)
{
	for(long i = 0; i < index.size(); ++i)
	{
		long pre_idx = index[i].start; long current_type = 1;
		long idx = index[i].start;
		while (idx <= index[i].end)
		{
			// percent_rank = ( rank() -1 ) / (total row -1 )
			if (!is_same(rankOptions.get_data()[pre_idx], rankOptions.get_data()[idx], measure_values))
			{
				pre_idx = idx;
				current_type = idx - index[i].start + 1;
			}
			double base = ( index[i].end - index[i].start == 0 ) ? 1 : index[i].end - index[i].start;
			copy_results(rankOptions, (current_type - 1.0 )/base, idx, get_data(result));
			idx++;
		};
	}
	return rankOptions.get_row();
};

uint32_t Util::mark_percrank_p(char *result, RankOptions& rankOptions, std::vector<Pair>& index, 
		std::vector<long>& measure_values, std::vector<double>& sum)
{
	tbb::parallel_for(size_t(0), index.size(), size_t(1) , [&](size_t i) {
			long pre_idx = index[i].start; long current_type = 1;
			long idx = index[i].start;
			while (idx <= index[i].end)
			{
			if (!is_same(rankOptions.get_data()[pre_idx], rankOptions.get_data()[idx], measure_values))
			{
			pre_idx = idx;
				current_type = idx - index[i].start + 1;
			}
			double base = ( index[i].end - index[i].start == 0 ) ? 1 : index[i].end - index[i].start;
			copy_results(rankOptions, (current_type - 1.0 )/base, idx, get_data(result));
			idx++;
			};
			});
	return rankOptions.get_row();
};

uint32_t Util::mark_denserank_p(char *result, RankOptions& rankOptions, std::vector<Pair>& index, 
		std::vector<long>& measure_values, std::vector<double>& sum)
{
	tbb::parallel_for(size_t(0), index.size(), size_t(1) , [&](size_t i) {
			long pre_idx = index[i].start; long current_type = 1;
			long idx = index[i].start;
			while (idx <= index[i].end)
			{
			if (!is_same(rankOptions.get_data()[pre_idx], rankOptions.get_data()[idx], measure_values))
			{
			pre_idx = idx;
			current_type++;
			}
			copy_results(rankOptions, current_type, idx, get_data(result));
			idx++;
			};
			});
	return rankOptions.get_row();
};


uint32_t Util::mark_denserank(char *result, RankOptions& rankOptions, std::vector<Pair>& index, 
		std::vector<long>& measure_values, std::vector<double>& sum)
{
	long part_idx = 0;
	for(long i = 0; i < index.size(); ++i)
	{
		long pre_idx = index[i].start; long current_type = 1;
		while (part_idx <= index[i].end)
		{
			if (!is_same(rankOptions.get_data()[pre_idx], 
						rankOptions.get_data()[part_idx], measure_values))
			{
				pre_idx = part_idx;
				current_type++;
			}
			copy_results(rankOptions, current_type, part_idx, get_data(result));
			part_idx++;
		};
	};
	return rankOptions.get_row();
};



uint32_t Util::mark_ntile(char *result, RankOptions& rankOptions, std::vector<Pair>& index, 
		std::vector<long>& measure_values, std::vector<double>& sum)
{
	long part_idx;
	for(long i = 0; i < index.size(); ++i)
	{
		long part = (index[i].end - index[i].start + 1) / rankOptions.range_vecs[0];
		long remains = (index[i].end - index[i].start + 1) % rankOptions.range_vecs[0];
		double current_value=0; long current_type = 0;
		part_idx = index[i].start;
		while (part_idx <= index[i].end)
		{
			current_value ++;
			copy_results(rankOptions, current_type, part_idx, get_data(result));
			if (current_value == part + (remains > 0))
			{
				current_value = 0;
				current_type++;
				remains--;
			};
			part_idx++;
		};
	};
	return part_idx;
};

uint32_t Util::mark_ntile_p(char *result, RankOptions& rankOptions, std::vector<Pair>& index, 
		std::vector<long>& measure_values, std::vector<double>& sum)
{
	tbb::parallel_for(size_t(0), index.size(), size_t(1) , [&](size_t i) {
			long part = (index[i].end - index[i].start + 1) / rankOptions.range_vecs[0];
			long remains = (index[i].end - index[i].start + 1) % rankOptions.range_vecs[0];
			double current_value=0; long current_type = 0;
			long part_idx = index[i].start;
			while (part_idx <= index[i].end)
			{
			current_value ++;
			copy_results(rankOptions, current_type, part_idx, get_data(result));
			if (current_value == part + (remains > 0))
			{
			current_value = 0;
			current_type++;
			remains--;
			};
			part_idx++;
			};
	});
	return rankOptions.get_row(); 
};


uint32_t Util::mark_tb_count(char *result, RankOptions& rankOptions, std::vector<Pair>& index, 
		std::vector<long>& measure_values, std::vector<double>& sum)
{
	long start = 0;
	long end = rankOptions.range_vecs[0];
	if ( rankOptions.range_vecs.size() == 2)
	{
		start = end - 1;
		end = rankOptions.range_vecs[1];
	};
	long result_idx = 0;
	//TODO: check end <= index[i] - part_idx + 1;
	for (long i = 0; i < index.size(); ++i )
	{
		if ( index[i].start + start > index[i].end ) 
		{
			continue;
		};
		int idx = index[i].start + start;
		long idx_end = index[i].start + end > index[i].end ? index[i].end + 1 : index[i].start + end;
		while ( idx < idx_end )
		{
			for (long j = 0; j < rankOptions.result_cols.size(); ++j)
			{
				get_data(result)[result_idx][j] = rankOptions.get_data()[idx][rankOptions.result_cols[j]];
			};
			result_idx++;
			idx++;
		}
		// with ties;
		while ( rankOptions.with_ties && idx <= index[i].end && 
				is_same(rankOptions.get_data()[idx], rankOptions.get_data()[idx - 1], measure_values))
		{
			for (long j = 0; j < rankOptions.result_cols.size(); ++j)
			{
				get_data(result)[result_idx][j] = rankOptions.get_data()[idx][rankOptions.result_cols[j]];
			};
			result_idx++;
			idx++;
		}
	};
	return (uint32_t) result_idx;
};


uint32_t Util::mark_tb_percent(char *result, RankOptions& rankOptions, std::vector<Pair>& index, 
		std::vector<long>& measure_values, std::vector<double>& sum)
{
	double start = 0;
	double end = rankOptions.range_vecs[0]/100.0;
	if ( rankOptions.range_vecs.size() == 2)
	{
		start = end;
		end = rankOptions.range_vecs[1]/100.0;
	};
	long result_idx = 0;
	for (long i = 0; i < index.size(); ++i )
	{
		long part_idx = index[i].start; 
		while ( part_idx <= index[i].end )
		{
			//include this one
			if ((part_idx - index[i].start +1)/(index[i].end - index[i].start + 1) >= start && 
					(part_idx - index[i].start + 1)/(index[i].end - index[i].start + 1) < end)
			{

				for (long j = 0; j < rankOptions.result_cols.size(); ++j)
				{
						get_data(result)[result_idx][j] = rankOptions.get_data()[part_idx][rankOptions.result_cols[j]];
				};
				result_idx++;
			}
			part_idx++;
			if ( (part_idx - index[i].start + 1)/(index[i].end - index[i].start + 1) > end)
			{
				while ( rankOptions.with_ties && part_idx <= index[i].end &&
						is_same(rankOptions.get_data()[part_idx], rankOptions.get_data()[part_idx - 1], measure_values))
				{
					for (long j = 0; j < rankOptions.result_cols.size(); ++j)
					{
						get_data(result)[result_idx][j] = rankOptions.get_data()[part_idx][rankOptions.result_cols[j]];
					};
					result_idx++;
					part_idx++;
				};
				break;
			};
		}
	};
	return (uint32_t) result_idx;
};

uint32_t Util::mark_tb_percent_on(char *result, RankOptions& rankOptions, std::vector<Pair>& index,
		std::vector<long>& measure_values, std::vector<double>& sum)
{
	double start = 0;
	double end = rankOptions.range_vecs[0]/100.0;
	if ( rankOptions.range_vecs.size() == 2)
	{
		start = end;
		end = rankOptions.range_vecs[1]/100.0;
	};
	uint32_t result_idx = 0;
	for (long i = 0; i < index.size(); ++i )
	{
		double running_total = 0;
		long part_idx = index[i].start; 
		while ( part_idx <= index[i].end )
		{
			running_total += rankOptions.get_data()[part_idx][rankOptions.percent_on_col];
			double percentage = sum[i] == 0.0 ?  0 : running_total/sum[i];
			if (percentage >= start && percentage < end)
			{
				for (long j = 0; j < rankOptions.result_cols.size(); ++j)
				{
					get_data(result)[result_idx][j] = rankOptions.get_data()[part_idx][rankOptions.result_cols[j]];
				};
				result_idx++;
			}
			part_idx++;
			if ( percentage > end)
			{
				while ( rankOptions.with_ties && part_idx <= index[i].end && 
						is_same(rankOptions.get_data()[part_idx], rankOptions.get_data()[part_idx - 1], measure_values))
				{
					for (long j = 0; j < rankOptions.result_cols.size(); ++j)
					{
						get_data(result)[result_idx][j] = rankOptions.get_data()[part_idx][rankOptions.result_cols[j]];
					};
					result_idx++;
					part_idx++;
				};
				break;
			};
		}
	};
	return result_idx;
};

void Util::psum(RankOptions& rankOptions, std::vector<Pair>& index)
{
	long idx = 1; 
	long start = 0;
	while (idx < rankOptions.get_row())
	{
		if (!is_same(rankOptions.get_data()[idx], rankOptions.get_data()[idx-1], rankOptions.dim_cols))
		{
			index.push_back(Pair(start, idx-1));
			start = idx;
		};
		idx++;
	};
	index.push_back(Pair(start,idx-1));
};

void Util::psum(RankOptions& rankOptions, std::vector<double>& sum, std::vector<Pair>& index)
{
	long idx = 1;  
	long start = 0;
	double tmp_sum = (rankOptions.get_data())[0][rankOptions.percent_on_col];
	while (idx < rankOptions.get_row())
	{
		if (!is_same(rankOptions.get_data()[idx], rankOptions.get_data()[idx-1], rankOptions.dim_cols))
		{   
			index.push_back(Pair(start, idx-1));
			sum.emplace_back(tmp_sum);
			tmp_sum = 0;
			start = idx;
		};  
		tmp_sum += rankOptions.get_data()[idx] [rankOptions.percent_on_col];
		idx++;
	}
	index.push_back(Pair(start,idx-1));sum.emplace_back(tmp_sum);
};

char* Util::new_new(const long row, const long measure_col_cnt, const long dim_col_cnt)
{
	//std::cout << "new_new: row " << row <<" ," << col << std::endl;
	//12 bytes padding for the header
	long col = measure_col_cnt + dim_col_cnt;
	char *storage = new char[ row * col * sizeof(double) + row * sizeof(double) + 12]; 
	uint32_t* header = (uint32_t*) storage;
	header[0] = (uint32_t) row;
	header[1] = (uint32_t) measure_col_cnt;
	header[2] = (uint32_t) dim_col_cnt;
	header += 3;
	double** data = (double **)header;
	for (int i = 0; i < row; ++i )
	{
		data[i] = (double*)(storage + 3 * sizeof(uint32_t) + row * sizeof(double) + i * col * sizeof(double)); 
	};  
	return storage;
};

char* Util::sort(RankOptions& rankOptions)
{
	long max_thread = 0;
	const long P_CUTOFF = 2000000;
	// convert measure_cols:INT_MIN
	std::vector<long> measures;
	for_each(rankOptions.measure_cols.begin(), rankOptions.measure_cols.end(), [&](long e)
			{
			if ( e == INT_MIN)
			measures.push_back(-0);
			else
			measures.emplace_back(e);		 

			});

	//for(auto i : measures)
	//	std::cout << "measures: " << i << std::endl;
	Functor cmp(rankOptions.dim_cols, measures);

	uint32_t *header = (uint32_t*)rankOptions.data;
	uint32_t row = header[0];
	uint32_t col_1 = header[1];
	uint32_t col_2 = header[2];
	//std::cout << "row is " << row << std::endl;
	header += 3;
	double ** data = (double **) header;
	// strip header pass data only
	//std::sort(data, data+row, cmp);
	//tbb::parallel_sort(data, data+row, [] (double * a, double * b) { return a[1] > b[1];});
	//tbb::parallel_sort(data, data+row, cmp);
  parallel_sample_sort(data, data+row, cmp);
	//parallel_merge_sort_tbb(data, data+row, cmp);
	std::vector<Pair> index; std::vector<double> sum; std::vector<long> measure_val;
	for(auto e : rankOptions.measure_cols)
	{
		if ( e == INT_MIN )
			measure_val.emplace_back(0);
		else
			measure_val.emplace_back( std::abs(e));
	};
	if (rankOptions.op == eTB_PERCENT_ON)
		psum(rankOptions, sum, index);
	else
		psum(rankOptions, index);

	// calculate rows
	if ( rankOptions.op == eTB_PERCENT)
	{
		if ( rankOptions.range_vecs.size() == 1)
			row = std::ceil(row * rankOptions.range_vecs[0]);  
		else
			row = std::ceil(row * rankOptions.range_vecs[1] - row * rankOptions.range_vecs[0]);  
	};
	char * result = new_new(row, rankOptions.measure_col_cnt, rankOptions.dim_col_cnt);
	uint32_t result_row = rankOptions.get_row() <  P_CUTOFF ? Util::op_map[rankOptions.op](result, rankOptions, index, measure_val, sum)
		: Util::op_map_parallel[rankOptions.op](result, rankOptions, index, measure_val, sum);
	set_row(result, result_row);

	delete[] rankOptions.data;
	return result;
};

/*********end of function/class **********/
constexpr Util::OP_FUNC Util::op_map_parallel[7];
constexpr Util::OP_FUNC Util::op_map[7];

char* PSort::sort(RankOptions& rankOptions)
{
	return Util::sort(rankOptions);
};

