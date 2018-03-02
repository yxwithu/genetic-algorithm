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
	void initPopulation();  //��ʼ����Ⱥ
	void crossover(); //����
	void mutation();  //����
	void estimationAndSelect(); //���ۡ�ѡ�񡢱���
	void run();  //����
	Individual& get_best_individual();  //�õ����Ÿ���
	~GA();

private:
	const unsigned size = 100;  //��Ⱥ��С
	const double cross_ratio = 0.7;  //�������
	const double sequence_mutation_ratio = 0.05;  //���б������
	const double separator_mutation_ratio = 0.05; //�ϵ�������
	const int iteration_num = 3000;  //��������
	vector<Individual> population;
	Individual best_individual;  //��Ⱥ�����ŵĸ���

	void crossover(Individual& individual1, Individual& individual2); //��������֮�佻��
	void sequence_random_swap(Individual& individual);  //���ѡȡSequence�������Ļ�����н��� 
	void separator_random_adjust(Individual& individual);  //���ѡȡseparator�������ϵ���е�����һ��һ����ֵͬ

	bool sort_by_score(const Individual& individual1, const Individual& individual2);  //���ڰ��÷�����
};

