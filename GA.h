#pragma once

#include "Individual.h"
#include <vector>
#include <algorithm>
#include <time.h>
#include <map>
#include <stdlib.h>
using namespace std;

class GA
{
public:
	GA();
	void initPopulation();  //初始化种群
	void crossover(); //交叉
	void mutation();  //变异
	void estimationAndSelect(); //评价、选择、保优
	void run();  //运行
	Individual& get_best_individual();  //得到最优个体
	~GA();

private:
	const unsigned size = 100;  //种群大小
	const double cross_ratio = 0.7;  //交叉概率
	const double sequence_mutation_ratio = 0.05;  //序列变异概率
	const double separator_mutation_ratio = 0.05; //断点变异概率
	const int iteration_num = 3000;  //迭代次数
	vector<Individual> population;
	Individual best_individual;  //种群中最优的个体

	void crossover(Individual& individual1, Individual& individual2); //两个个体之间交叉
	void sequence_random_swap(Individual& individual);  //随机选取Sequence中两个的基因进行交换 
	void separator_random_adjust(Individual& individual);  //随机选取separator中两个断点进行调整，一加一减相同值

	bool sort_by_score(const Individual& individual1, const Individual& individual2);  //用于按得分排序
};

