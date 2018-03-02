#include "GA.h"


GA::GA()
{
}

Individual& GA::get_best_individual()
{
	return best_individual;
}

void GA::run()
{
	initPopulation();  //初始化
	estimationAndSelect();  //评价选择保优

	for (int i = 0; i < iteration_num; ++i){
		crossover();  //交叉
		mutation();  //变异
		estimationAndSelect();  //评价选择保优
	}
}

void GA::initPopulation()
{
	if (population.size() > 0){
		population.clear();   //多次运行先清除
	}
	//初始化size个个体，加入到种群中，不需要计算得分
	for (int i = 0; i < size; ++i){
		Individual individual;
		individual.Init_Idividuals();
		population.push_back(individual);
	}
	Individual temp;
	temp.Init_Idividuals();
	best_individual = temp;//初始化最优个体也是个随机个体
}

void GA::estimationAndSelect()
{
	//评价
	vector<double> fitness;  //适应度，存放累计适应度
	double fitness_sum = 0;  //用作归一化
	int max_idx = 0;  //适应度最强的个体索引
	double max_value = 0;
	for (int i = 0; i < size; ++i){
		population[i].getObj();
		double value = 1.0 / population[i].ObjValue;
		fitness.push_back(value);
		fitness_sum += value;

		if (value > max_value){
			max_idx = i;
			max_value = value;
		}
	}
	best_individual = population[max_idx]; //用于保优
	
	//为轮盘赌选择做准备，适应度转换成被选择的概率
	double sum_ratio = 0;
	for (int i = 0; i < size; ++i){
		double ratio = fitness[i] / fitness_sum;
		sum_ratio += ratio;
		fitness[i] = sum_ratio;
	}

	//根据适应度进行选择，二分查找，同时找到适应度最小的个体用于替换成保优个体
	vector<Individual> select_res;
	srand((unsigned)time(NULL));

	int min_idx = 0;  //新种群中适应度最低的个体
	double min_value = 1;

	for (int i = 0; i < size; ++i){
		double ratio = rand() / double(RAND_MAX);
		int left = 0, right = size - 1, mid = 0;
		while (left < right){
			mid = (left + right) / 2;
			//mid对应的概率区间
			double lowValue = mid == 0 ? 0 : fitness[mid - 1];
			double highValue = fitness[mid];
			if (ratio >= lowValue && ratio < highValue){
				break;
			}
			else if (ratio < lowValue){
				right = mid - 1;
			}
			else{
				left = mid + 1;
			}
		}
		if (left >= right){
			mid = left;
		}
		select_res.push_back(population[mid]);

		double value = mid == 0 ? fitness[mid] : fitness[mid] - fitness[mid - 1];
		if (value < min_value){
			min_idx = i;
			min_value = value;
		}
	}

	select_res[min_idx] = best_individual;  //保优

	//更新种群
	population.clear();
	for (int i = 0; i < size; ++i){
		population.push_back(select_res[i]);
	}
}

bool GA::sort_by_score(const Individual& individual1, const Individual& individual2)
{
	return individual1.ObjValue < individual2.ObjValue;
}


void GA::crossover()
{
	//随机分组，两两交叉
	int idx[100];
	for (int i = 0; i < 100; ++i){
		idx[i] = i;
	}
	//打乱数组，再按序取两两一组，就等于是随机分组了
	random_shuffle(&idx[0], &idx[99]);

	//开始进行交叉
	for (int i = 0; i < 50; ++i){
		int idx1 = idx[2 * i];
		int idx2 = idx[2 * i + 1];
		crossover(population.at(idx1), population.at(idx2));
	}
}

void GA::crossover(Individual& individual1, Individual& individual2)
{
	srand((unsigned)time(NULL));
	if (rand() / double(RAND_MAX) > cross_ratio) return;  //不发生交叉
	
	//对序列基因数组
	int len = individual1.SequenceChromosome.size();
	vector<int> zero_marks;  //为0的部分，为1的部分不需要变动，只需要变动为0的部分
	for (int i = 0; i < len; ++i){
		if (rand() % 2 == 0){
			zero_marks.push_back(i);
		}
	}

	map<int, int> g1_idx_map, g2_idx_map;  //得到父代序列基因值与位置的映射关系，value1-index1, value2-index2
	for (int i = 0; i < len; ++i){
		g1_idx_map[individual1.SequenceChromosome[i]] = i;
		g2_idx_map[individual2.SequenceChromosome[i]] = i;
	}

	map<int, int> f1_value_map, f2_value_map;  //index2-value1, index1-value2
	for (int mark : zero_marks){
		int g1_seq_value = individual1.SequenceChromosome[mark]; //g1中0标记位置上的值
		int g2_seq_index = g2_idx_map[g1_seq_value];  //这个值在g2中的位置
		f1_value_map[g2_seq_index] = g1_seq_value;  //利用map将g1中的值按g2中的顺序排序

		int g2_seq_value = individual2.SequenceChromosome[mark];
		int g1_seq_index = g1_idx_map[g2_seq_value];
		f2_value_map[g1_seq_index] = g2_seq_value;
	}

	map<int, int>::iterator iter;  //将染色体标记为0的部分按排好的顺序赋值
	int i = 0;
	for (iter = f1_value_map.begin(); iter != f1_value_map.end(); i++,iter++){
		individual1.SequenceChromosome[zero_marks[i]] = iter->second;
	}
	for (i = 0, iter = f2_value_map.begin(); iter != f2_value_map.end(); i++, iter++){
		individual2.SequenceChromosome[zero_marks[i]] = iter->second; 
	}

	//对断点基因数组，进行交换
	if (rand() / double(RAND_MAX) <= 0.5){
		individual1.SeparatorChromosome.swap(individual2.SeparatorChromosome);
	}
}

void GA::mutation(){
	//对种群中的每一个个体进行变异
	srand((unsigned)time(NULL));
	
	//先对序列数组进行变异
	for (int i = 0; i < size; ++i){ 
		if (rand() / double(RAND_MAX) <= sequence_mutation_ratio){
			sequence_random_swap(population[i]);  //序列数组随机交换两个位置上的值
		}
	}

	//再对断点序列进行变异
	for (int i = 0; i < size; ++i){
		if (rand() / double(RAND_MAX) <= separator_mutation_ratio){
			separator_random_adjust(population[i]);  //随机改变断点数组两个位置上的值，保证总数不变
		}
	}
}

void GA::sequence_random_swap(Individual& individual)
{
	srand((unsigned)time(NULL));
	int idx1, idx2;
	int len = individual.SequenceChromosome.size();
	idx1 = rand() % len;
	idx2 = rand() % len;
	while (idx1 == idx2)//随机生成两个不同索引
	{
		idx2 = rand() % len;
	}

	int temp;
	temp = individual.SequenceChromosome[idx1];
	individual.SequenceChromosome[idx1] = individual.SequenceChromosome[idx2];
	individual.SequenceChromosome[idx2] = temp;
}

void GA::separator_random_adjust(Individual& individual)
{
	srand((unsigned)time(NULL));
	int len = individual.SeparatorChromosome.size();
	if (len == 1) return;  //一行的无法调整
	int idx1, idx2;  //两个变化的index
	idx1 = rand() % len;  //减
	while (individual.SeparatorChromosome[idx1] == 0){
		idx1 = rand() % len;
	}
	idx2 = rand() % len;  //加
	while (idx1 == idx2){
		idx2 = rand() % len;
	}
	//改变的值
	int value = rand() % individual.SeparatorChromosome[idx1];
	individual.SeparatorChromosome[idx1] -= value;
	individual.SeparatorChromosome[idx2] += value;
}

GA::~GA()
{
}
