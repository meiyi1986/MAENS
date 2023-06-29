
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "functions.h"

#if INF != 2100000000
#define INF 2100000000
#define DEPOT 1
#define DUMMY_CYCLE 0
#define MAX_TASK_NUM 380
#define MAX_TASK_TAG_LENGTH 381
#define MAX_NODE_TAG_LENGTH 141
#define MAX_TASK_ROUTE_LENGTH 193
#define MAX_NODE_ROUTE_LENGTH 300
#define MAX_TASK_SEQ_LENGTH 250
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void mod_dijkstra()
{
	int i, j, k, m, minimum;
	
	for (i = 1; i <= vertex_num; i++)
	{
		for (j = 1; j <= vertex_num; j++)
		{
			if (j == i)
				continue;
			
			shortest_path[i][j][0] = 1;
			shortest_path[i][j][1] = i;
			min_cost[i][j] = INF;
		}
	}

	int mark[MAX_NODE_TAG_LENGTH], dist[MAX_NODE_TAG_LENGTH], dist1[MAX_NODE_TAG_LENGTH], nearest_neighbor[MAX_NODE_TAG_LENGTH];

	for (i = 1; i <= vertex_num; i++)
	{
		mark[i] = 1;
		
		for (j = 1; j <= vertex_num; j++)
		{
			if (j == i)
				continue;
			
			mark[j] = 0;
			dist[j] = trav_cost[i][j];
			dist1[j] = dist[j];
		}
		
		for (k = 1; k < vertex_num; k++)
		{
			minimum = INF;
			nearest_neighbor[0] = 0;
			
			for (j = 1; j <= vertex_num; j++)
			{
				if (mark[j])
					continue;
					
				if (dist1[j] == INF)
					continue;
					
				if (dist1[j] < minimum)
					minimum = dist1[j];
			}
			
			if (minimum == INF)
				continue;
			
			for (j = 1; j <= vertex_num; j++)
			{
				if (mark[j])
					continue;
					
				if (dist1[j] == minimum)
				{
					nearest_neighbor[0] ++;
					nearest_neighbor[nearest_neighbor[0]] = j;
				}
			}
			
			int v = nearest_neighbor[1];
			dist1[v] = INF;
			mark[v] = 1;
			
			if (shortest_path[i][v][0] == 0 || (shortest_path[i][v][0] > 0 && shortest_path[i][v][shortest_path[i][v][0]] != v))
			{
				shortest_path[i][v][0] ++;
				shortest_path[i][v][shortest_path[i][v][0]] = v;
			}
				
			for (j = 1; j <= vertex_num; j++)
			{
				if (mark[j])
					continue;
					
				if (minimum+trav_cost[v][j] < dist[j])
				{
					dist[j] = minimum+trav_cost[v][j];
					dist1[j] = minimum+trav_cost[v][j];
					for (m = 0; m <= shortest_path[i][v][0]; m++)
					{
						shortest_path[i][j][m] = shortest_path[i][v][m];
					}
				}
			}
			
			for (j = 1; j <= vertex_num; j++)
			{
				if (j == i)
					continue;
				
				min_cost[i][j] = dist[j];
			}
		}		
	}
	
	for (i = 1; i <= vertex_num; i++)
	{
		for (j = 1; j <= vertex_num; j++)
		{
			if (shortest_path[i][j][0] == 1)
				shortest_path[i][j][0] = 0;
		}
	}
	
	for (i = 1; i <= vertex_num; i++)
	{
		shortest_path[i][i][0] = 1;
		shortest_path[i][i][1] = i;
		min_cost[i][i] = 0;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void path_scanning(individual *ps_indi, const task *inst_tasks, int *serve_mark)
{
	int i, k;
	int serve_task_num = 0;
	for (i = req_edge_num+1; i <= task_num; i++)
	{
		if (serve_mark[i])
			serve_task_num ++;
	}
	
	int load, trial, mindist;

	int unserved_task[MAX_TASK_TAG_LENGTH], candi_task[MAX_TASK_TAG_LENGTH], nearest_task[MAX_TASK_TAG_LENGTH];
	int nearest_isol_task[MAX_TASK_TAG_LENGTH], nearest_inci_task[MAX_TASK_TAG_LENGTH], sel_task[MAX_TASK_TAG_LENGTH];
	int current_task, next_task;

	int positions[MAX_TASK_SEQ_LENGTH];
	
	individual tmp_indi1, tmp_indi2, tmp_indi3, tmp_indi4, tmp_indi5;

	ps_indi->total_cost = INF;

	int dep_dist[MAX_TASK_TAG_LENGTH];
	double yield[MAX_TASK_TAG_LENGTH];
	int max_dep_dist, min_dep_dist;
	double max_yield, min_yield;
	
	for (i = 1; i <= task_num; i++)
	{
		if (!serve_mark[i])
			continue;
		
		dep_dist[i] = min_cost[inst_tasks[i].tail_node][DEPOT];
		yield[i] = 1.0*inst_tasks[i].demand/inst_tasks[i].serv_cost;
	}
	
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    tmp_indi1.sequence[0] = 1;
	tmp_indi1.sequence[1] = 0;
	tmp_indi1.route_seg_load[0] = 0;
	
	unserved_task[0] = 0;
	for(i = 1; i <= task_num; i++)
	{
		if (!serve_mark[i])
			continue;
		
		unserved_task[0] ++;
		unserved_task[unserved_task[0]] = i;
	}
	
	load = 0;
	trial = 0;
	while (trial < serve_task_num)
	{
		current_task = tmp_indi1.sequence[tmp_indi1.sequence[0]];
		
		candi_task[0] = 0;
		
		for (i = 1; i <= unserved_task[0]; i++)
		{
			if (inst_tasks[unserved_task[i]].demand <= capacity-load)
			{
				candi_task[0] ++;
				candi_task[candi_task[0]] = unserved_task[i];
			}
		}
		
		if (candi_task[0] == 0)
		{
  			tmp_indi1.sequence[0] ++;
	  		tmp_indi1.sequence[tmp_indi1.sequence[0]] = 0;
		  	tmp_indi1.route_seg_load[0] ++;
			tmp_indi1.route_seg_load[tmp_indi1.route_seg_load[0]] = load;
			load = 0;
			continue;
		}
		
		mindist = INF;
		nearest_task[0] = 0;
		
		for (i = 1; i <= candi_task[0]; i++)
		{
			if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] < mindist)
			{
				mindist = min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node];
				nearest_task[0] = 1;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
			else if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] == mindist)
			{
				nearest_task[0] ++;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
		}
		
		nearest_isol_task[0] = 0;
		nearest_inci_task[0] = 0;
		
		for (i = 1; i <= nearest_task[0]; i++)
		{
			if (inst_tasks[nearest_task[i]].tail_node == 1)
			{
				nearest_inci_task[0] ++;
				nearest_inci_task[nearest_inci_task[0]] = nearest_task[i];
			}
			else
			{
				nearest_isol_task[0] ++;
				nearest_isol_task[nearest_isol_task[0]] = nearest_task[i];
			}
		}
		
		if (nearest_isol_task[0] == 0)
		{
			memcpy(nearest_isol_task, nearest_inci_task, (nearest_inci_task[0]+1)*sizeof(int));
		}
		
		max_dep_dist = -1;
		sel_task[0] = 0;
		for (i = 1; i <= nearest_isol_task[0]; i++)
		{
			if (dep_dist[nearest_isol_task[i]] > max_dep_dist)
			{
				max_dep_dist = dep_dist[nearest_isol_task[i]];
				sel_task[0] = 1;
				sel_task[sel_task[0]] = nearest_isol_task[i];
			}
			else if (dep_dist[nearest_isol_task[i]] == max_dep_dist)
			{
				sel_task[0] ++;
				sel_task[sel_task[0]] = nearest_isol_task[i];
			}
		}
		
		//k = rand_choose(sel_task[0]);
		k = 1;
		next_task = sel_task[k];
		
		trial ++;
		tmp_indi1.sequence[0] ++;
		tmp_indi1.sequence[tmp_indi1.sequence[0]] = next_task;
		load += inst_tasks[next_task].demand;
		find_ele_positions(positions, unserved_task, next_task);
		delete_element(unserved_task, positions[1]);
		
		if (inst_tasks[next_task].inverse > 0)
		{
			find_ele_positions(positions, unserved_task, inst_tasks[next_task].inverse);
			delete_element(unserved_task, positions[1]);
		}
	}
	
	tmp_indi1.sequence[0] ++;
	tmp_indi1.sequence[tmp_indi1.sequence[0]] = 0;
	tmp_indi1.route_seg_load[0] ++;
	tmp_indi1.route_seg_load[tmp_indi1.route_seg_load[0]] = load;
		
	tmp_indi1.total_cost = get_task_seq_total_cost(tmp_indi1.sequence, inst_tasks);
	tmp_indi1.total_vio_load = get_total_vio_load(tmp_indi1.route_seg_load);
	
	if (tmp_indi1.total_cost < ps_indi->total_cost)
	{
		memcpy(ps_indi->sequence, tmp_indi1.sequence, (tmp_indi1.sequence[0]+1)*sizeof(int));
		memcpy(ps_indi->route_seg_load, tmp_indi1.route_seg_load, (tmp_indi1.route_seg_load[0]+1)*sizeof(int));
		ps_indi->total_cost = tmp_indi1.total_cost;
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	
	tmp_indi2.sequence[0] = 1;
	tmp_indi2.sequence[1] = 0;
	tmp_indi2.route_seg_load[0] = 0;
	
	unserved_task[0] = 0;
	for(i = 1; i <= task_num; i++)
	{
		if (!serve_mark[i])
			continue;
		
		unserved_task[0] ++;
		unserved_task[unserved_task[0]] = i;
	}
	
	load = 0;
	trial = 0;
	while (trial < serve_task_num)
	{
		current_task = tmp_indi2.sequence[tmp_indi2.sequence[0]];
		
		candi_task[0] = 0;
		
		for (i = 1; i <= unserved_task[0]; i++)
		{
			if (inst_tasks[unserved_task[i]].demand <= capacity-load)
			{
				candi_task[0] ++;
				candi_task[candi_task[0]] = unserved_task[i];
			}
		}
		
		if (candi_task[0] == 0)
		{
			tmp_indi2.sequence[0] ++;
	  		tmp_indi2.sequence[tmp_indi2.sequence[0]] = 0;
		  	tmp_indi2.route_seg_load[0] ++;
			tmp_indi2.route_seg_load[tmp_indi2.route_seg_load[0]] = load;
			load = 0;
			continue;
		}
		
		mindist = INF;
		nearest_task[0] = 0;
		
		for (i = 1; i <= candi_task[0]; i++)
		{
			if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] < mindist)
			{
				mindist = min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node];
				nearest_task[0] = 1;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
			else if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] == mindist)
			{
				nearest_task[0] ++;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
		}
		
		nearest_isol_task[0] = 0;
		nearest_inci_task[0] = 0;
		
		for (i = 1; i <= nearest_task[0]; i++)
		{
			if (inst_tasks[nearest_task[i]].tail_node == 1)
			{
				nearest_inci_task[0] ++;
				nearest_inci_task[nearest_inci_task[0]] = nearest_task[i];
			}
			else
			{
				nearest_isol_task[0] ++;
				nearest_isol_task[nearest_isol_task[0]] = nearest_task[i];
			}
		}
		
		if (nearest_isol_task[0] == 0)
		{
			memcpy(nearest_isol_task, nearest_inci_task, (nearest_inci_task[0]+1)*sizeof(int));
		}
		
		min_dep_dist = INF;
		sel_task[0] = 0;
		for (i = 1; i <= nearest_isol_task[0]; i++)
		{
			if (dep_dist[nearest_isol_task[i]] < min_dep_dist)
			{
				min_dep_dist = dep_dist[nearest_isol_task[i]];
				sel_task[0] = 1;
				sel_task[sel_task[0]] = nearest_isol_task[i];
			}
			else if (dep_dist[nearest_isol_task[i]] == min_dep_dist)
			{
				sel_task[0] ++;
				sel_task[sel_task[0]] = nearest_isol_task[i];
			}
		}
		
		//k = rand_choose(sel_task[0]);
		k = 1;
		next_task = sel_task[k];
		
		trial ++;
		tmp_indi2.sequence[0] ++;
		tmp_indi2.sequence[tmp_indi2.sequence[0]] = next_task;
		load += inst_tasks[next_task].demand;
		find_ele_positions(positions, unserved_task, next_task);
		delete_element(unserved_task, positions[1]);
		
		if (inst_tasks[next_task].inverse > 0)
		{
			find_ele_positions(positions, unserved_task, inst_tasks[next_task].inverse);
			delete_element(unserved_task, positions[1]);
		}
	}
	
	tmp_indi2.sequence[0] ++;
	tmp_indi2.sequence[tmp_indi2.sequence[0]] = 0;
	tmp_indi2.route_seg_load[0] ++;
	tmp_indi2.route_seg_load[tmp_indi2.route_seg_load[0]] = load;
		
	tmp_indi2.total_cost = get_task_seq_total_cost(tmp_indi2.sequence, inst_tasks);
	tmp_indi2.total_vio_load = get_total_vio_load(tmp_indi2.route_seg_load);
	
	if (tmp_indi2.total_cost < ps_indi->total_cost)
	{
		memcpy(ps_indi->sequence, tmp_indi2.sequence, (tmp_indi2.sequence[0]+1)*sizeof(int));
		memcpy(ps_indi->route_seg_load, tmp_indi2.route_seg_load, (tmp_indi2.route_seg_load[0]+1)*sizeof(int));
		ps_indi->total_cost = tmp_indi2.total_cost;
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	
	tmp_indi3.sequence[0] = 1;
	tmp_indi3.sequence[1] = 0;
	tmp_indi3.route_seg_load[0] = 0;
	
	unserved_task[0] = 0;
	for(i = 1; i <= task_num; i++)
	{
		if (!serve_mark[i])
			continue;
		
		unserved_task[0] ++;
		unserved_task[unserved_task[0]] = i;
	}
	
	load = 0;
	trial = 0;
	while (trial < serve_task_num)
	{
		current_task = tmp_indi3.sequence[tmp_indi3.sequence[0]];
		
		candi_task[0] = 0;
		
		for (i = 1; i <= unserved_task[0]; i++)
		{
			if (inst_tasks[unserved_task[i]].demand <= capacity-load)
			{
				candi_task[0] ++;
				candi_task[candi_task[0]] = unserved_task[i];
			}
		}
		
		if (candi_task[0] == 0)
		{
  			tmp_indi3.sequence[0] ++;
	  		tmp_indi3.sequence[tmp_indi3.sequence[0]] = 0;
		  	tmp_indi3.route_seg_load[0] ++;
			tmp_indi3.route_seg_load[tmp_indi3.route_seg_load[0]] = load;
			load = 0;
			continue;
		}
		
		mindist = INF;
		nearest_task[0] = 0;
		
		for (i = 1; i <= candi_task[0]; i++)
		{
			if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] < mindist)
			{
				mindist = min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node];
				nearest_task[0] = 1;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
			else if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] == mindist)
			{
				nearest_task[0] ++;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
		}
		
		nearest_isol_task[0] = 0;
		nearest_inci_task[0] = 0;
		
		for (i = 1; i <= nearest_task[0]; i++)
		{
			if (inst_tasks[nearest_task[i]].tail_node == 1)
			{
				nearest_inci_task[0] ++;
				nearest_inci_task[nearest_inci_task[0]] = nearest_task[i];
			}
			else
			{
				nearest_isol_task[0] ++;
				nearest_isol_task[nearest_isol_task[0]] = nearest_task[i];
			}
		}
		
		if (nearest_isol_task[0] == 0)
		{
			memcpy(nearest_isol_task, nearest_inci_task, (nearest_inci_task[0]+1)*sizeof(int));
		}
		
		max_yield = -1;
		sel_task[0] = 0;
		for (i = 1; i <= nearest_isol_task[0]; i++)
		{
			if (yield[nearest_isol_task[i]] > max_yield)
			{
				max_yield = yield[nearest_isol_task[i]];
				sel_task[0] = 1;
				sel_task[sel_task[0]] = nearest_isol_task[i];
			}
			else if (yield[nearest_isol_task[i]] == max_yield)
			{
				sel_task[0] ++;
				sel_task[sel_task[0]] = nearest_isol_task[i];
			}
		}
		
		//k = rand_choose(sel_task[0]);
		k = 1;
		next_task = sel_task[k];
		
		trial ++;
		tmp_indi3.sequence[0] ++;
		tmp_indi3.sequence[tmp_indi3.sequence[0]] = next_task;
		load += inst_tasks[next_task].demand;
		find_ele_positions(positions, unserved_task, next_task);
		delete_element(unserved_task, positions[1]);
		
		if (inst_tasks[next_task].inverse > 0)
		{
			find_ele_positions(positions, unserved_task, inst_tasks[next_task].inverse);
			delete_element(unserved_task, positions[1]);
		}
	}
	
	tmp_indi3.sequence[0] ++;
	tmp_indi3.sequence[tmp_indi3.sequence[0]] = 0;
	tmp_indi3.route_seg_load[0] ++;
	tmp_indi3.route_seg_load[tmp_indi3.route_seg_load[0]] = load;
		
	tmp_indi3.total_cost = get_task_seq_total_cost(tmp_indi3.sequence, inst_tasks);
	tmp_indi3.total_vio_load = get_total_vio_load(tmp_indi3.route_seg_load);
	
	if (tmp_indi3.total_cost < ps_indi->total_cost)
	{
		memcpy(ps_indi->sequence, tmp_indi3.sequence, (tmp_indi3.sequence[0]+1)*sizeof(int));
		memcpy(ps_indi->route_seg_load, tmp_indi3.route_seg_load, (tmp_indi3.route_seg_load[0]+1)*sizeof(int));
		ps_indi->total_cost = tmp_indi3.total_cost;
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	
	tmp_indi4.sequence[0] = 1;
	tmp_indi4.sequence[1] = 0;
	tmp_indi4.route_seg_load[0] = 0;
	
	unserved_task[0] = 0;
	for(i = 1; i <= task_num; i++)
	{
		if (!serve_mark[i])
			continue;
		
		unserved_task[0] ++;
		unserved_task[unserved_task[0]] = i;
	}
	
	load = 0;
	trial = 0;
	while (trial < serve_task_num)
	{
		current_task = tmp_indi4.sequence[tmp_indi4.sequence[0]];
		
		candi_task[0] = 0;
		
		for (i = 1; i <= unserved_task[0]; i++)
		{
			if (inst_tasks[unserved_task[i]].demand <= capacity-load)
			{
				candi_task[0] ++;
				candi_task[candi_task[0]] = unserved_task[i];
			}
		}
		
		if (candi_task[0] == 0)
		{
  			tmp_indi4.sequence[0] ++;
	  		tmp_indi4.sequence[tmp_indi4.sequence[0]] = 0;
		  	tmp_indi4.route_seg_load[0] ++;
			tmp_indi4.route_seg_load[tmp_indi4.route_seg_load[0]] = load;
			load = 0;
			continue;
		}
		
		mindist = INF;
		nearest_task[0] = 0;
		
		for (i = 1; i <= candi_task[0]; i++)
		{
			if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] < mindist)
			{
				mindist = min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node];
				nearest_task[0] = 1;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
			else if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] == mindist)
			{
				nearest_task[0] ++;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
		}
		
		nearest_isol_task[0] = 0;
		nearest_inci_task[0] = 0;
		
		for (i = 1; i <= nearest_task[0]; i++)
		{
			if (inst_tasks[nearest_task[i]].tail_node == 1)
			{
				nearest_inci_task[0] ++;
				nearest_inci_task[nearest_inci_task[0]] = nearest_task[i];
			}
			else
			{
				nearest_isol_task[0] ++;
				nearest_isol_task[nearest_isol_task[0]] = nearest_task[i];
			}
		}
		
		if (nearest_isol_task[0] == 0)
		{
			memcpy(nearest_isol_task, nearest_inci_task, (nearest_inci_task[0]+1)*sizeof(int));
		}
		
		min_yield = INF;
		sel_task[0] = 0;
		for (i = 1; i <= nearest_isol_task[0]; i++)
		{
			if (yield[nearest_isol_task[i]] < min_yield)
			{
				min_yield = yield[nearest_isol_task[i]];
				sel_task[0] = 1;
				sel_task[sel_task[0]] = nearest_isol_task[i];
			}
			else if (yield[nearest_isol_task[i]] == min_yield)
			{
				sel_task[0] ++;
				sel_task[sel_task[0]] = nearest_isol_task[i];
			}
		}
		
		//k = rand_choose(sel_task[0]);
		k = 1;
		next_task = sel_task[k];
		
		trial ++;
		tmp_indi4.sequence[0] ++;
		tmp_indi4.sequence[tmp_indi4.sequence[0]] = next_task;
		load += inst_tasks[next_task].demand;
		find_ele_positions(positions, unserved_task, next_task);
		delete_element(unserved_task, positions[1]);
		
		if (inst_tasks[next_task].inverse > 0)
		{
			find_ele_positions(positions, unserved_task, inst_tasks[next_task].inverse);
			delete_element(unserved_task, positions[1]);
		}
	}
	
	tmp_indi4.sequence[0] ++;
	tmp_indi4.sequence[tmp_indi4.sequence[0]] = 0;
	tmp_indi4.route_seg_load[0] ++;
	tmp_indi4.route_seg_load[tmp_indi4.route_seg_load[0]] = load;
		
	tmp_indi4.total_cost = get_task_seq_total_cost(tmp_indi4.sequence, inst_tasks);
	tmp_indi4.total_vio_load = get_total_vio_load(tmp_indi4.route_seg_load);
	
	if (tmp_indi4.total_cost < ps_indi->total_cost)
	{
		memcpy(ps_indi->sequence, tmp_indi4.sequence, (tmp_indi4.sequence[0]+1)*sizeof(int));
		memcpy(ps_indi->route_seg_load, tmp_indi4.route_seg_load, (tmp_indi4.route_seg_load[0]+1)*sizeof(int));
		ps_indi->total_cost = tmp_indi4.total_cost;
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
    tmp_indi5.sequence[0] = 1;
	tmp_indi5.sequence[1] = 0;
	tmp_indi5.route_seg_load[0] = 0;
	
	unserved_task[0] = 0;
	for(i = 1; i <= task_num; i++)
	{
		if (!serve_mark[i])
			continue;
		
		unserved_task[0] ++;
		unserved_task[unserved_task[0]] = i;
	}
	
	load = 0;
	trial = 0;
	while (trial < serve_task_num)
	{
		current_task = tmp_indi5.sequence[tmp_indi5.sequence[0]];
		
		candi_task[0] = 0;
		
		for (i = 1; i <= unserved_task[0]; i++)
		{
			if (inst_tasks[unserved_task[i]].demand <= capacity-load)
			{
				candi_task[0] ++;
				candi_task[candi_task[0]] = unserved_task[i];
			}
		}
		
		if (candi_task[0] == 0)
		{
  			tmp_indi5.sequence[0] ++;
	  		tmp_indi5.sequence[tmp_indi5.sequence[0]] = 0;
		  	tmp_indi5.route_seg_load[0] ++;
			tmp_indi5.route_seg_load[tmp_indi5.route_seg_load[0]] = load;
			load = 0;
			continue;
		}
		
		mindist = INF;
		nearest_task[0] = 0;
		
		for (i = 1; i <= candi_task[0]; i++)
		{
			if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] < mindist)
			{
				mindist = min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node];
				nearest_task[0] = 1;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
			else if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] == mindist)
			{
				nearest_task[0] ++;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
		}
		
		nearest_isol_task[0] = 0;
		nearest_inci_task[0] = 0;
		
		for (i = 1; i <= nearest_task[0]; i++)
		{
			if (inst_tasks[nearest_task[i]].tail_node == 1)
			{
				nearest_inci_task[0] ++;
				nearest_inci_task[nearest_inci_task[0]] = nearest_task[i];
			}
			else
			{
				nearest_isol_task[0] ++;
				nearest_isol_task[nearest_isol_task[0]] = nearest_task[i];
			}
		}
		
		if (nearest_isol_task[0] == 0)
		{
			memcpy(nearest_isol_task, nearest_inci_task, (nearest_inci_task[0]+1)*sizeof(int));
		}
		
		if (load < capacity/2)
		{
			max_dep_dist = -1;
			sel_task[0] = 0;
			for (i = 1; i <= nearest_isol_task[0]; i++)
			{
				if (dep_dist[nearest_isol_task[i]] > max_dep_dist)
				{
					max_dep_dist = dep_dist[nearest_isol_task[i]];
					sel_task[0] = 1;
					sel_task[sel_task[0]] = nearest_isol_task[i];
				}
				else if (dep_dist[nearest_isol_task[i]] == max_dep_dist)
				{
					sel_task[0] ++;
					sel_task[sel_task[0]] = nearest_isol_task[i];
				}
			}
		}
		else
		{
			min_dep_dist = INF;
			sel_task[0] = 0;
			for (i = 1; i <= nearest_isol_task[0]; i++)
			{
				if (dep_dist[nearest_isol_task[i]] < min_dep_dist)
				{
					min_dep_dist = dep_dist[nearest_isol_task[i]];
					sel_task[0] = 1;
					sel_task[sel_task[0]] = nearest_isol_task[i];
				}
				else if (dep_dist[nearest_isol_task[i]] == min_dep_dist)
				{
					sel_task[0] ++;
					sel_task[sel_task[0]] = nearest_isol_task[i];
				}
			}
		}
		
		//k = rand_choose(sel_task[0]);
		k = 1;
		next_task = sel_task[k];
		
		trial ++;
		tmp_indi5.sequence[0] ++;
		tmp_indi5.sequence[tmp_indi5.sequence[0]] = next_task;
		load += inst_tasks[next_task].demand;
		find_ele_positions(positions, unserved_task, next_task);
		delete_element(unserved_task, positions[1]);
		
		if (inst_tasks[next_task].inverse > 0)
		{
			find_ele_positions(positions, unserved_task, inst_tasks[next_task].inverse);
			delete_element(unserved_task, positions[1]);
		}
	}
	
	tmp_indi5.sequence[0] ++;
	tmp_indi5.sequence[tmp_indi5.sequence[0]] = 0;
	tmp_indi5.route_seg_load[0] ++;
	tmp_indi5.route_seg_load[tmp_indi5.route_seg_load[0]] = load;
	
	tmp_indi5.total_cost = get_task_seq_total_cost(tmp_indi5.sequence, inst_tasks);
	tmp_indi5.total_vio_load = get_total_vio_load(tmp_indi5.route_seg_load);
	
	if (tmp_indi5.total_cost < ps_indi->total_cost)
	{
		memcpy(ps_indi->sequence, tmp_indi5.sequence, (tmp_indi5.sequence[0]+1)*sizeof(int));
		memcpy(ps_indi->route_seg_load, tmp_indi5.route_seg_load, (tmp_indi5.route_seg_load[0]+1)*sizeof(int));
		ps_indi->total_cost = tmp_indi5.total_cost;
	}

	ps_indi->total_vio_load = 0;
	//get_route_seg_length(ps_indi->route_seg_length, ps_indi->sequence, inst_tasks);
	//ps_indi->max_length = max(ps_indi->route_seg_length);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void rand_scanning(individual *rs_indi, const task *inst_tasks, int *serve_mark)
{
	int i, k;
	int serve_task_num = 0;
	for (i = req_edge_num+1; i <= task_num; i++)
	{
		if (serve_mark[i])
			serve_task_num ++;
	}
	
	int load, trial, mindist;

	int unserved_task[MAX_TASK_TAG_LENGTH], candi_task[MAX_TASK_TAG_LENGTH], nearest_task[MAX_TASK_TAG_LENGTH];
	int current_task, next_task;

	int positions[MAX_TASK_SEQ_LENGTH];
	
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    rs_indi->sequence[0] = 1;
	rs_indi->sequence[1] = 0;
	rs_indi->route_seg_load[0] = 0;
	
	unserved_task[0] = 0;
	for(i = 1; i <= task_num; i++)
	{
		if (!serve_mark[i])
			continue;
		
		unserved_task[0] ++;
		unserved_task[unserved_task[0]] = i;
	}
	
	load = 0;
	trial = 0;
	while (trial < serve_task_num)
	{
		current_task = rs_indi->sequence[rs_indi->sequence[0]];
		
		candi_task[0] = 0;
		
		for (i = 1; i <= unserved_task[0]; i++)
		{
			if (inst_tasks[unserved_task[i]].demand <= capacity-load)
			{
				candi_task[0] ++;
				candi_task[candi_task[0]] = unserved_task[i];
			}
		}
		
		if (candi_task[0] == 0)
		{
  			rs_indi->sequence[0] ++;
	  		rs_indi->sequence[rs_indi->sequence[0]] = 0;
		  	rs_indi->route_seg_load[0] ++;
			rs_indi->route_seg_load[rs_indi->route_seg_load[0]] = load;
			load = 0;
			continue;
		}

		mindist = INF;
		nearest_task[0] = 0;
		
		for (i = 1; i <= candi_task[0]; i++)
		{
			if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] < mindist)
			{
				mindist = min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node];
				nearest_task[0] = 1;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
			else if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] == mindist)
			{
				nearest_task[0] ++;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
		}
		
		k = rand_choose(nearest_task[0]);
		next_task = nearest_task[k];
		
		//k = rand_choose(candi_task[0]);
		//next_task = candi_task[k];
		
		trial ++;
		rs_indi->sequence[0] ++;
		rs_indi->sequence[rs_indi->sequence[0]] = next_task;
		load += inst_tasks[next_task].demand;
		find_ele_positions(positions, unserved_task, next_task);
		delete_element(unserved_task, positions[1]);
		
		if (inst_tasks[next_task].inverse > 0)
		{
			find_ele_positions(positions, unserved_task, inst_tasks[next_task].inverse);
			delete_element(unserved_task, positions[1]);
		}
	}
	
	rs_indi->sequence[0] ++;
	rs_indi->sequence[rs_indi->sequence[0]] = 0;
	rs_indi->route_seg_load[0] ++;
	rs_indi->route_seg_load[rs_indi->route_seg_load[0]] = load;
		
	rs_indi->total_cost = get_task_seq_total_cost(rs_indi->sequence, inst_tasks);
	rs_indi->total_vio_load = get_total_vio_load(rs_indi->route_seg_load);
	//get_route_seg_length(rs_indi->route_seg_length, rs_indi->sequence, inst_tasks);
	//rs_indi->max_length = max(rs_indi->route_seg_length);
}