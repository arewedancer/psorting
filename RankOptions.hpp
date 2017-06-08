#ifndef RANKOPTIONS_H
#define RANKOPTIONS_H
#include <vector>
#define Rank_Result_Index 50001
typedef enum
{eTB_COUNT=0,eTB_PERCENT,eTB_PERCENT_ON, eRANK,eDENSERANK,ePERCRANK,eNTILE} SORT_OP;

struct RankOptions
{
	std::vector<long>& dim_cols;
	std::vector<long>& measure_cols;
	std::vector<unsigned long>& result_cols;
	int64_t percent_on_col;
	bool with_ties;
	std::vector<long>& range_vecs;
	SORT_OP op;
	char* data; // input data pointer;
	int32_t measure_col_cnt;  // final result set measure column count; must be put into buffer header
	int32_t dim_col_cnt; // final result set dimension column count; must be put into buffer header

	explicit RankOptions(std::vector<long>& dim_cols_, std::vector<long>& measure_cols_, std::vector<unsigned long>& result_cols_,
			int64_t percent_on_col_, bool with_ties_, std::vector<long>& range_vecs_, SORT_OP op_, char* data_,
			int32_t measure_col_cnt_, int32_t dim_col_cnt_):dim_cols(dim_cols_), measure_cols(measure_cols_), result_cols(result_cols_),
			percent_on_col(percent_on_col_), with_ties(with_ties_), range_vecs(range_vecs_), op(op_), 
			data(data_), measure_col_cnt(measure_col_cnt_), dim_col_cnt(dim_col_cnt_){};

	inline double ** get_data()
	{
		char * return_data = data + 3 * sizeof(uint32_t);
		return (double **) return_data;
	};
	inline uint32_t get_row()
	{
		return ((uint32_t*)data)[0];
	};
	inline uint32_t get_measure_col()
	{
		return ((uint32_t*)data)[1];
	}
	inline uint32_t get_dim_col()
	{
		return ((uint32_t*)data)[2];
	}
};
struct PSort
{
	static char* sort(RankOptions& RankOptions);
};

#endif
