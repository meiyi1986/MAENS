#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "functions.h"

int vertex_num;
int req_edge_num;
int req_arc_num;
int nonreq_edge_num;
int nonreq_arc_num;
int task_num;
int total_arc_num;
int vehicle_num;
int capacity;

int DEPOT;

int trav_cost[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
int serve_cost[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
int shortest_path[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
int min_cost[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];

int popsize = 30; // population size can be set here

char dummy_string[50];

char input_file[100] = "instance/example.dat"; // the instance can be changed here
int lower_bound = 3548; // the lower bound of the instance can be set here

task inst_tasks[MAX_TASK_TAG_LENGTH];
arc inst_arcs[MAX_ARCS_TAG_LENGTH];

individual pop[MAX_TOTALSIZE];
individual best_fsb_solution;

int main(void)
{
	/* input and preprocessing */

	FILE *fp;

	req_arc_num = 0;
	nonreq_arc_num = 0;

	fp = fopen(input_file, "r");

	while (1)
	{
		fscanf(fp, "%s", dummy_string); // get the instance name
		//printf("dummy_string = %s\n", dummy_string);

		if (strcmp(dummy_string, "VERTICES") == 0)
		{
			fscanf(fp, "%s", dummy_string);
			fscanf(fp, "%d", &vertex_num);
		}
		else if (strcmp(dummy_string, "ARISTAS_REQ") == 0)
		{
			fscanf(fp, "%s", dummy_string);
			fscanf(fp, "%d", &req_edge_num);
		}
		else if (strcmp(dummy_string, "ARISTAS_NOREQ") == 0)
		{
			fscanf(fp, "%s", dummy_string);
			fscanf(fp, "%d", &nonreq_edge_num);
		}
		else if (strcmp(dummy_string, "VEHICULOS") == 0)
		{
			fscanf(fp, "%s", dummy_string);
			fscanf(fp, "%d", &vehicle_num);
		}
		else if (strcmp(dummy_string, "CAPACIDAD") == 0)
		{
			fscanf(fp, "%s", dummy_string);
			fscanf(fp, "%d", &capacity);
		}
		else if (strcmp(dummy_string, "LISTA_ARISTAS_REQ") == 0)
		{
			task_num = 2*req_edge_num+req_arc_num;
			total_arc_num = task_num+2*nonreq_edge_num+nonreq_arc_num;

			fscanf(fp, "%s", dummy_string);

			for (int i = 1; i <= req_edge_num; i++)
			{
				fscanf(fp, "%s", dummy_string); // get the left parenthesis "("
				fscanf(fp, "%d,", &inst_tasks[i].head_node);
				fscanf(fp, "%d)", &inst_tasks[i].tail_node);
				fscanf(fp, "%s", dummy_string); // get the right parenthesis ")" and the string "coste"
				fscanf(fp, "%d", &inst_tasks[i].serv_cost);
				fscanf(fp, "%s", dummy_string); // get the string "demanda"
				fscanf(fp, "%d", &inst_tasks[i].demand);
				inst_tasks[i].dead_cost = inst_tasks[i].serv_cost;
				inst_tasks[i].inverse = i+req_edge_num;

				inst_tasks[i+req_edge_num].head_node = inst_tasks[i].tail_node;
				inst_tasks[i+req_edge_num].tail_node = inst_tasks[i].head_node;
				inst_tasks[i+req_edge_num].dead_cost = inst_tasks[i].dead_cost;
				inst_tasks[i+req_edge_num].serv_cost = inst_tasks[i].serv_cost;
				inst_tasks[i+req_edge_num].demand = inst_tasks[i].demand;
				inst_tasks[i+req_edge_num].inverse = i;

				inst_arcs[i].head_node = inst_tasks[i].head_node;
				inst_arcs[i].tail_node = inst_tasks[i].tail_node;
				inst_arcs[i].trav_cost = inst_tasks[i].dead_cost;
				inst_arcs[i+req_edge_num].head_node = inst_arcs[i].tail_node;
				inst_arcs[i+req_edge_num].tail_node = inst_arcs[i].head_node;
				inst_arcs[i+req_edge_num].trav_cost = inst_arcs[i].trav_cost;
			}
		}
		else if (strcmp(dummy_string, "LISTA_ARISTAS_NOREQ") == 0)
		{
			fscanf(fp, "%s", dummy_string);

			for (int i = task_num+1; i <= task_num+nonreq_edge_num; i++)
			{
				fscanf(fp, "%s", dummy_string); // get the left parenthesis "("
				fscanf(fp, "%d,", &inst_arcs[i].head_node);
				fscanf(fp, "%d)", &inst_arcs[i].tail_node);
				fscanf(fp, "%s", dummy_string); // get the right parenthesis ")" and the string "coste"
				fscanf(fp, "%d", &inst_arcs[i].trav_cost);

				inst_arcs[i+nonreq_edge_num].head_node = inst_arcs[i].tail_node;
				inst_arcs[i+nonreq_edge_num].tail_node = inst_arcs[i].head_node;
				inst_arcs[i+nonreq_edge_num].trav_cost = inst_arcs[i].trav_cost;
			}
		}
		else if (strcmp(dummy_string, "DEPOSITO") == 0)
		{
			fscanf(fp, "%s", dummy_string);
			fscanf(fp, "%d", &DEPOT);
			break;
		}
	}

	fclose(fp);

	inst_tasks[DUMMY_CYCLE].tail_node = DEPOT;
	inst_tasks[DUMMY_CYCLE].head_node = DEPOT;
	inst_tasks[DUMMY_CYCLE].dead_cost = 0;
	inst_tasks[DUMMY_CYCLE].serv_cost = 0;
	inst_tasks[DUMMY_CYCLE].demand = 0;
	inst_tasks[DUMMY_CYCLE].inverse = DUMMY_CYCLE;
	inst_arcs[DUMMY_CYCLE].tail_node = DEPOT;
	inst_arcs[DUMMY_CYCLE].head_node = DEPOT;
	inst_arcs[DUMMY_CYCLE].trav_cost = 0;

	for (int i = 1; i <= vertex_num; i++)
	{
		for (int j = 1; j <= vertex_num; j++)
		{
			trav_cost[i][j] = INF;
			serve_cost[i][j] = 0;
		}
	}

	for (int i = 1; i <= total_arc_num; i++)
	{
		trav_cost[inst_arcs[i].head_node][inst_arcs[i].tail_node] = inst_arcs[i].trav_cost;
	}

	for (int i = 1; i <= task_num; i++)
	{
		serve_cost[inst_tasks[i].head_node][inst_tasks[i].tail_node] = inst_tasks[i].serv_cost;
	}

	mod_dijkstra();

	/* get output_file */

	char output_file[100] = "output.dat"; // the name of the output file can be set here

	/* random seed */

	clock_t start = clock();

	int tm;
	tm = time(NULL);
	//tm = 1239002157;
	srand(tm);

	//printf("tm = %d\n", tm);

	fp = fopen(output_file, "w");

	fprintf(fp, "The random seed is %d.\n", tm);
	fprintf(fp, "There are %d vertices, %d tasks, and the capacities of vehicles is %d.\n\n", vertex_num, req_edge_num+req_arc_num, capacity);

	fclose(fp);

	/* initialization */

	best_fsb_solution.total_cost = INF;

	int used;

	int tmp_popsize = 0;
	while (tmp_popsize < popsize)
	{
		int trial = 0;
		individual init_indi;
		while (trial < M_trial)
		{
			trial ++;
			int serve_mark[MAX_TASK_TAG_LENGTH];
			memset(serve_mark, 0, sizeof(serve_mark));
			serve_mark[0] = task_num;
			for (int i = 1; i <= task_num; i++)
			{
				serve_mark[i] = 1;
			}

			//path_scanning(&init_indi, inst_tasks, serve_mark);
			rand_scanning(&init_indi, inst_tasks, serve_mark);

			used = 0;
			for (int i = 0; i < tmp_popsize; i++)
			{
				if (init_indi.total_cost == pop[i].total_cost && init_indi.total_vio_load == pop[i].total_vio_load)
				{
					used = 1;
					break;
				}
			}

			if (!used)
				break;
		}

		if (trial == M_trial && used)
			break;

		indi_copy(&pop[tmp_popsize], &init_indi);
		tmp_popsize ++;			

		if (init_indi.total_vio_load == 0 && init_indi.total_cost < best_fsb_solution.total_cost)
		{
			indi_copy(&best_fsb_solution, &init_indi);
		}
	}
	popsize = tmp_popsize;
	/*
	for (int i = 0; i < popsize; i++)
	{
	get_predsucc(pop[i].pred, pop[i].succ, pop[i].sequence, inst_tasks);
	}

	for (int i = 0; i < popsize; i++)
	{
	for (int j = 0; j < popsize; j++)
	{
	int dist = calc_distance(&pop[i], &pop[j], inst_tasks);
	printf("%d  ", dist);
	}
	printf("\n");
	}
	*/
	/* searching phase */

	int ite = 0, wite = 0;
	individual parent1, parent2, xed_child, mted_child, child;

	int offsize = 6*popsize; // the number of the generated offsprings
	int totalsize = popsize+offsize; // the number of the sum of the current population and the offsprings

	while (ite < M_ite)
	{
		//for (int i = 0; i < popsize; i++)
		//{
		//	printf("%d ", pop[i].total_cost);
		//}
		//printf("\n");
		ite ++;
		wite ++;

		//if (ite%25 == 0)
		//{
		//	printf("ite = %d, best cost = %d\n", ite, best_fsb_solution.total_cost);
		//}

		/* generate offsprings */

		int ptr = popsize;
		while (ptr < totalsize)
		{
			child.total_cost = 0; // child does not exist

			int par_id1, par_id2;
			//tour_selection(par_id1, par_id2, pop); // parent selection
			rand_selection(par_id1, par_id2, pop);
			indi_copy(&parent1, &pop[par_id1]);
			indi_copy(&parent2, &pop[par_id2]);

			SBX(&xed_child, &parent1, &parent2, inst_tasks); // crossover

			if (xed_child.total_vio_load == 0 && xed_child.total_cost < best_fsb_solution.total_cost)
			{
				indi_copy(&best_fsb_solution, &xed_child);
				wite = 0;

				fp = fopen(output_file, "a");

				fprintf(fp, "crossover of ite = %d\n", ite);
				fprintf(fp, "sequence\n");
				for (int i = 1; i <= best_fsb_solution.sequence[0]; i++)
				{
					fprintf(fp, "%d  ", best_fsb_solution.sequence[i]);
				}
				fprintf(fp, "\ntotal_cost = %d\n\n", best_fsb_solution.total_cost);

				fclose(fp);
			}

			used = 0;
			for (int i = 0; i < ptr; i++)
			{
				if (i == par_id1 || i == par_id2)
					continue;

				if (xed_child.total_cost == pop[i].total_cost && xed_child.total_vio_load == pop[i].total_vio_load)
				{
					used = 1;
					break;
				}
			}

			if (!used)
			{
				indi_copy(&child, &xed_child);
			}

			/* mutation with probability M_PROB */

			double random = 1.0*rand()/RAND_MAX;
			if (random < M_PROB)
			{
				/*
				int tmp_tc = best_fsb_solution.total_cost;
				*/
				lns_mut(&mted_child, &xed_child, inst_tasks);

				//if (mted_child.total_vio_load == 0 && mted_child.total_cost < best_fsb_solution.total_cost)
				//{
				//	indi_copy(&best_fsb_solution, &mted_child);
				//	printf("new best cost = %d\n", best_fsb_solution.total_cost);
				//}
				/*
				if (best_fsb_solution.total_cost < tmp_tc)
				{
				wite = 0;

				fp = fopen(output_file, "a");

				fprintf(fp, "local search of ite = %d\n", ite);
				fprintf(fp, "sequence\n");
				for (int i = 1; i <= best_fsb_solution.sequence[0]; i++)
				{
				fprintf(fp, "%d  ", best_fsb_solution.sequence[i]);
				}
				fprintf(fp, "\ntotal_cost = %d\n\n", best_fsb_solution.total_cost);

				fclose(fp);
				}
				*/
				used = 0;
				for (int i = 0; i < ptr; i++)
				{
					if (i == par_id1 || i == par_id2)
						continue;

					if (mted_child.total_cost == pop[i].total_cost && mted_child.total_vio_load == pop[i].total_vio_load)
					{
						used = 1;
						break;
					}
				}

				if (!used)
				{
					indi_copy(&child, &mted_child);
				}
			}

			if (child.total_cost == parent1.total_cost && child.total_vio_load == parent1.total_vio_load)
			{
				indi_copy(&pop[par_id1], &child);
			}
			else if (child.total_cost == parent2.total_cost && child.total_vio_load == parent2.total_vio_load)
			{
				indi_copy(&pop[par_id2], &child);
			}
			else if (child.total_cost > 0)
			{
				indi_copy(&pop[ptr], &child);
				ptr ++;
			}

			if (best_fsb_solution.total_cost == lower_bound)
				break;
		}

		if (best_fsb_solution.total_cost == lower_bound)
			break;

		// stochastic ranking

		double Pf = 0.45;
		individual tmp_indi;

		for (int i = 0; i < totalsize; i++)
		{
			for (int j = 0; j < i; j++)
			{
				double random = 1.0*rand()/RAND_MAX;
				if ((pop[j].total_vio_load == 0 && pop[j+1].total_vio_load == 0) || random < Pf)
				{
					if (pop[j].total_cost > pop[j+1].total_cost)
					{
						indi_copy(&tmp_indi, &pop[j]);
						indi_copy(&pop[j], &pop[j+1]);
						indi_copy(&pop[j+1], &tmp_indi);
					}
				}
				else
				{
					if (pop[j].total_vio_load > pop[j+1].total_vio_load)
					{
						indi_copy(&tmp_indi, &pop[j]);
						indi_copy(&pop[j], &pop[j+1]);
						indi_copy(&pop[j+1], &tmp_indi);
					}
				}
			}
		}
	}

	/*for (int i = 0; i < popsize; i++)
	{
	get_predsucc(pop[i].pred, pop[i].succ, pop[i].sequence, inst_tasks);
	}

	for (int i = 0; i < popsize; i++)
	{
	for (int j = 0; j < popsize; j++)
	{
	int dist = calc_distance(&pop[i], &pop[j], inst_tasks);
	printf("%d  ", dist);
	}
	printf("\n");
	}*/

	clock_t finish = clock();
	double duration = (double)(finish-start)/CLOCKS_PER_SEC;

	fp = fopen(output_file, "a");

	fprintf(fp, "The best sequence is:\n");
	for (int i = 1; i <= best_fsb_solution.sequence[0]; i++)
	{
		fprintf(fp, "%d ", best_fsb_solution.sequence[i]);
	}
	fprintf(fp, "\nThe best total cost is %d.\n", best_fsb_solution.total_cost);
	fprintf(fp, "The duration time is %f seconds.\n", duration);
	fprintf(fp, "The number of iterations is %d.\n", ite);

	fclose(fp);
}