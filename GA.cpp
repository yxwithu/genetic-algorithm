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
	initPopulation();  //��ʼ��
	estimationAndSelect();  //����ѡ����

	for (int i = 0; i < iteration_num; ++i){
		crossover();  //����
		mutation();  //����
		estimationAndSelect();  //����ѡ����
	}
}

void GA::initPopulation()
{
	if (population.size() > 0){
		population.clear();   //������������
	}
	//��ʼ��size�����壬���뵽��Ⱥ�У�����Ҫ����÷�
	for (int i = 0; i < size; ++i){
		Individual individual;
		individual.Init_Idividuals();
		population.push_back(individual);
	}
	Individual temp;
	temp.Init_Idividuals();
	best_individual = temp;//��ʼ�����Ÿ���Ҳ�Ǹ��������
}

void GA::estimationAndSelect()
{
	//����
	vector<double> fitness;  //��Ӧ�ȣ�����ۼ���Ӧ��
	double fitness_sum = 0;  //������һ��
	int max_idx = 0;  //��Ӧ����ǿ�ĸ�������
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
	best_individual = population[max_idx]; //���ڱ���
	
	//Ϊ���̶�ѡ����׼������Ӧ��ת���ɱ�ѡ��ĸ���
	double sum_ratio = 0;
	for (int i = 0; i < size; ++i){
		double ratio = fitness[i] / fitness_sum;
		sum_ratio += ratio;
		fitness[i] = sum_ratio;
	}

	//������Ӧ�Ƚ���ѡ�񣬶��ֲ��ң�ͬʱ�ҵ���Ӧ����С�ĸ��������滻�ɱ��Ÿ���
	vector<Individual> select_res;
	srand((unsigned)time(NULL));

	int min_idx = 0;  //����Ⱥ����Ӧ����͵ĸ���
	double min_value = 1;

	for (int i = 0; i < size; ++i){
		double ratio = rand() / double(RAND_MAX);
		int left = 0, right = size - 1, mid = 0;
		while (left < right){
			mid = (left + right) / 2;
			//mid��Ӧ�ĸ�������
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

	select_res[min_idx] = best_individual;  //����

	//������Ⱥ
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
	//������飬��������
	int idx[100];
	for (int i = 0; i < 100; ++i){
		idx[i] = i;
	}
	//�������飬�ٰ���ȡ����һ�飬�͵��������������
	random_shuffle(&idx[0], &idx[99]);

	//��ʼ���н���
	for (int i = 0; i < 50; ++i){
		int idx1 = idx[2 * i];
		int idx2 = idx[2 * i + 1];
		crossover(population.at(idx1), population.at(idx2));
	}
}

void GA::crossover(Individual& individual1, Individual& individual2)
{
	srand((unsigned)time(NULL));
	if (rand() / double(RAND_MAX) > cross_ratio) return;  //����������
	
	//�����л�������
	int len = individual1.SequenceChromosome.size();
	vector<int> zero_marks;  //Ϊ0�Ĳ��֣�Ϊ1�Ĳ��ֲ���Ҫ�䶯��ֻ��Ҫ�䶯Ϊ0�Ĳ���
	for (int i = 0; i < len; ++i){
		if (rand() % 2 == 0){
			zero_marks.push_back(i);
		}
	}

	map<int, int> g1_idx_map, g2_idx_map;  //�õ��������л���ֵ��λ�õ�ӳ���ϵ��value1-index1, value2-index2
	for (int i = 0; i < len; ++i){
		g1_idx_map[individual1.SequenceChromosome[i]] = i;
		g2_idx_map[individual2.SequenceChromosome[i]] = i;
	}

	map<int, int> f1_value_map, f2_value_map;  //index2-value1, index1-value2
	for (int mark : zero_marks){
		int g1_seq_value = individual1.SequenceChromosome[mark]; //g1��0���λ���ϵ�ֵ
		int g2_seq_index = g2_idx_map[g1_seq_value];  //���ֵ��g2�е�λ��
		f1_value_map[g2_seq_index] = g1_seq_value;  //����map��g1�е�ֵ��g2�е�˳������

		int g2_seq_value = individual2.SequenceChromosome[mark];
		int g1_seq_index = g1_idx_map[g2_seq_value];
		f2_value_map[g1_seq_index] = g2_seq_value;
	}

	map<int, int>::iterator iter;  //��Ⱦɫ����Ϊ0�Ĳ��ְ��źõ�˳��ֵ
	int i = 0;
	for (iter = f1_value_map.begin(); iter != f1_value_map.end(); i++,iter++){
		individual1.SequenceChromosome[zero_marks[i]] = iter->second;
	}
	for (i = 0, iter = f2_value_map.begin(); iter != f2_value_map.end(); i++, iter++){
		individual2.SequenceChromosome[zero_marks[i]] = iter->second; 
	}

	//�Զϵ�������飬���н���
	if (rand() / double(RAND_MAX) <= 0.5){
		individual1.SeparatorChromosome.swap(individual2.SeparatorChromosome);
	}
}

void GA::mutation(){
	//����Ⱥ�е�ÿһ��������б���
	srand((unsigned)time(NULL));
	
	//�ȶ�����������б���
	for (int i = 0; i < size; ++i){ 
		if (rand() / double(RAND_MAX) <= sequence_mutation_ratio){
			sequence_random_swap(population[i]);  //�������������������λ���ϵ�ֵ
		}
	}

	//�ٶԶϵ����н��б���
	for (int i = 0; i < size; ++i){
		if (rand() / double(RAND_MAX) <= separator_mutation_ratio){
			separator_random_adjust(population[i]);  //����ı�ϵ���������λ���ϵ�ֵ����֤��������
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
	while (idx1 == idx2)//�������������ͬ����
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
	if (len == 1) return;  //һ�е��޷�����
	int idx1, idx2;  //�����仯��index
	idx1 = rand() % len;  //��
	while (individual.SeparatorChromosome[idx1] == 0){
		idx1 = rand() % len;
	}
	idx2 = rand() % len;  //��
	while (idx1 == idx2){
		idx2 = rand() % len;
	}
	//�ı��ֵ
	int value = rand() % individual.SeparatorChromosome[idx1];
	individual.SeparatorChromosome[idx1] -= value;
	individual.SeparatorChromosome[idx2] += value;
}

GA::~GA()
{
}
