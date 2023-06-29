
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "functions.h"

int task_routes[MAX_SEG_TAG_LENGTH][MAX_TASK_SEG_LENGTH];

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void rand_selection(int &id1, int &id2, individual *pop)
/* pop has been sorted increasingly already */
{
	int k1, k2;
	int candi[MAX_POPSIZE+1];
	candi[0] = popsize;
	for (int i = 1; i <= popsize; i++)
	{
		candi[i] = i-1;
	}

	k1 = rand_choose(candi[0]);
	id1 = candi[k1];
	delete_element(candi, k1);
	k2 = rand_choose(candi[0]);
	id2 = candi[k2];
	//printf("id1 = %d, id2 = %d, popsize = %d\n", id1, id2, popsize);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void tour_selection(int &id1, int &id2, individual *pop)
/* pop has been sorted increasingly already */
{
	int k1, k2;
	int candi1[MAX_POPSIZE+1], candi2[MAX_POPSIZE+1];
	candi1[0] = popsize;
	for (int i = 1; i <= popsize; i++)
	{
		candi1[i] = i-1;
	}
	memcpy(candi2, candi1, sizeof(candi1));

	k1 = rand_choose(candi1[0]);
	delete_element(candi1, k1);
	k2 = rand_choose(candi1[0]);
	if (k1 < k2)
	{
		id1 = candi1[k1];
		delete_element(candi2, k1);
	}
	else
	{
		id1 = candi1[k2];
		delete_element(candi2, k2);
	}

	k1 = rand_choose(candi2[0]);
	delete_element(candi2, k1);
	k2 = rand_choose(candi2[0]);
	if (k1 < k2)
	{
		id2 = candi2[k1];
	}
	else
	{
		id2 = candi2[k2];
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void SBX(individual *xed_child, individual *p1, individual *p2, const task *inst_tasks)
{
	int task_seq1[MAX_TASK_SEQ_LENGTH], task_seq2[MAX_TASK_SEQ_LENGTH];
	int task_segs1[MAX_SEG_TAG_LENGTH][MAX_TASK_SEQ_LENGTH], task_segs2[MAX_SEG_TAG_LENGTH][MAX_TASK_SEQ_LENGTH];

	int positions[MAX_SEG_TAG_LENGTH];

	memcpy(task_seq1, p1->sequence, sizeof(p1->sequence));
	memcpy(task_seq2, p2->sequence, sizeof(p2->sequence));
	find_ele_positions(positions, task_seq1, 0);

	int i, j, k, l, m, n, w, k1, k2, NO_LeftTasks, IVLoad;
	int S1[MAX_TASK_SEQ_LENGTH], S2[MAX_TASK_SEQ_LENGTH], SubPath1[MAX_TASK_SEQ_LENGTH], SubPath2[MAX_TASK_SEQ_LENGTH];
	int Routes1[MAX_SEG_TAG_LENGTH][MAX_TASK_SEQ_LENGTH], Routes2[MAX_SEG_TAG_LENGTH][MAX_TASK_SEQ_LENGTH];
	int XCLds[MAX_SEG_TAG_LENGTH], LeftTasks[MAX_TASK_SEQ_LENGTH], Positions[MAX_SEG_TAG_LENGTH];
	
	struct CandSelection
	{
		int RouteID;
		int Pos;
	};
	
	struct CandSelection CandSLCTList1[MAX_TASK_SEQ_LENGTH], CandSLCTList2[MAX_TASK_SEQ_LENGTH];
	int LLength1 = 0, LLength2 = 0;
	
	memcpy(S1, p1->sequence, sizeof(p1->sequence));
	memcpy(S2, p2->sequence, sizeof(p2->sequence));
		
	find_ele_positions(Positions, S1, 0);
	Routes1[0][0] = Positions[0]-1;
	for (i = 1; i < Positions[0]; i++)
	{
		copy_sub_array(Routes1[i], S1, Positions[i], Positions[i+1]);
	}
	
	find_ele_positions(Positions, S2, 0);
	Routes2[0][0] = Positions[0]-1;
	for (i = 1; i < Positions[0]; i++)
	{
		copy_sub_array(Routes2[i], S2, Positions[i], Positions[i+1]);
	}
	
	memcpy(XCLds, p1->route_seg_load, sizeof(p1->route_seg_load));
	
	for (i = 1; i <= Routes1[0][0]; i++)
	{
		for (j = 2; j < Routes1[i][0]; j++)
		{
			LLength1 ++;
			CandSLCTList1[LLength1].RouteID = i;
			CandSLCTList1[LLength1].Pos = j;
		}
	}
	
	for (i = 1; i <= Routes2[0][0]; i++)
	{
		for (j = 2; j < Routes2[i][0]; j++)
		{
			LLength2 ++;
			CandSLCTList2[LLength2].RouteID = i;
			CandSLCTList2[LLength2].Pos = j;
		}
	}
	
	k1 = rand_choose(LLength1);
	k2 = rand_choose(LLength2);
	
	copy_sub_array(SubPath1, Routes1[CandSLCTList1[k1].RouteID], 1, CandSLCTList1[k1].Pos);
	copy_sub_array(SubPath2, Routes2[CandSLCTList2[k2].RouteID], CandSLCTList2[k2].Pos, Routes2[CandSLCTList2[k2].RouteID][0]);	
	copy_sub_array(LeftTasks, Routes1[CandSLCTList1[k1].RouteID], CandSLCTList1[k1].Pos+1, Routes1[CandSLCTList1[k1].RouteID][0]-1);
	
	int Checked[MAX_TASK_SEQ_LENGTH];
	memset(Checked, 0, sizeof(Checked));
	
	for (i = 1; i < SubPath2[0]; i++)
	{
		if (Checked[i])
			continue;
		
		for (j = SubPath1[0]; j > 1; j--)
		{
			if (SubPath1[j] == SubPath2[i] || SubPath1[j] == inst_tasks[SubPath2[i]].inverse)
			{
				delete_element(SubPath1, j);
				Checked[i] = 1;
				break;
			}
		}
	}
	
	for (i = 1; i < SubPath2[0]; i++)
	{
		if (Checked[i])
			continue;
		
		for (j = LeftTasks[0]; j > 0; j--)
		{
			if (LeftTasks[j] == SubPath2[i] || LeftTasks[j] == inst_tasks[SubPath2[i]].inverse)
			{
				delete_element(LeftTasks, j);
				Checked[i] = 1;
				break;
			}
		}
	}
	
	for (i = 1; i < SubPath2[0]; i++)
	{
		if (Checked[i])
			continue;
		
		for (j = 1; j <= Routes1[0][0]; j++)
		{
			if (j == CandSLCTList1[k1].RouteID)
			//if (j == ChangedRouteID)
				continue;
			
			for (k = Routes1[j][0]; k > 1; k--)
			{
				if (Routes1[j][k] == SubPath2[i] || Routes1[j][k] == inst_tasks[SubPath2[i]].inverse)
				{
					delete_element(Routes1[j], k);
					XCLds[j] -= inst_tasks[SubPath2[i]].demand;
					Checked[i] = 1;
					break;
				}
			}
			
			if (Checked[i])
				break;
		}
	}
	
	link_array(SubPath1, SubPath2);
	memcpy(Routes1[CandSLCTList1[k1].RouteID], SubPath1, sizeof(SubPath1));
	XCLds[CandSLCTList1[k1].RouteID] = 0;
	for (i = 2; i < Routes1[CandSLCTList1[k1].RouteID][0]; i++)
	{
		XCLds[CandSLCTList1[k1].RouteID] += inst_tasks[Routes1[CandSLCTList1[k1].RouteID][i]].demand;
	}
	
	NO_LeftTasks = LeftTasks[0];
	
	// insert left tasks
	
	struct Insert
	{
		int InsertedTask;
		int InsertRouteID;
		int InsertPos;
		int InsertCost;
		int InsertVioLoad;
	};
	
	struct Insert CandInsertions[6000];
	int NO_CandInsertions;
	
	struct Insert ParetoSetInsertions[6000];
	int ParetoSetSize;
	int Out[6000], Add;
	
	struct Insert BestInsertion;

	int NO_CurrRoutes = 0;
	for (int j = 1; j <= Routes1[0][0]; j++)
	{
		if (Routes1[j][0] == 2)
			continue;

		NO_CurrRoutes ++;
	}
	
	for (n = 1; n <= NO_LeftTasks; n++)
	{
		NO_CandInsertions = 0;
		ParetoSetSize = 0;
	  
		for (j = 1; j <= Routes1[0][0]; j++)
		{
			if (Routes1[j][0] == 2)
				continue;

			if (XCLds[j] > capacity)
			{
				IVLoad = inst_tasks[LeftTasks[n]].demand;
			}
			else if (XCLds[j] > capacity-inst_tasks[LeftTasks[n]].demand)
			{
				IVLoad = XCLds[j]+inst_tasks[LeftTasks[n]].demand-capacity;
			}
			else 
			{
				IVLoad = 0;
			}
   		
			for (k = 2; k <= Routes1[j][0]; k++)
			{
				NO_CandInsertions ++;
				CandInsertions[NO_CandInsertions].InsertedTask = LeftTasks[n];
				CandInsertions[NO_CandInsertions].InsertRouteID = j;
				CandInsertions[NO_CandInsertions].InsertPos = k;
				CandInsertions[NO_CandInsertions].InsertCost = min_cost[inst_tasks[Routes1[j][k-1]].tail_node][inst_tasks[LeftTasks[n]].head_node]+
					min_cost[inst_tasks[LeftTasks[n]].tail_node][inst_tasks[Routes1[j][k]].head_node]-
					min_cost[inst_tasks[Routes1[j][k-1]].tail_node][inst_tasks[Routes1[j][k]].head_node];
				CandInsertions[NO_CandInsertions].InsertVioLoad = IVLoad;
   			
				Out[0] = 0;
				Add = 1;
   			
				for (m = 1; m <= ParetoSetSize; m++)
				{
					if (CandInsertions[NO_CandInsertions].InsertCost > ParetoSetInsertions[m].InsertCost
						&& CandInsertions[NO_CandInsertions].InsertVioLoad > ParetoSetInsertions[m].InsertVioLoad)
					{
						Add = 0;
						break;
					}
					else if (CandInsertions[NO_CandInsertions].InsertCost < ParetoSetInsertions[m].InsertCost
						&& CandInsertions[NO_CandInsertions].InsertVioLoad < ParetoSetInsertions[m].InsertVioLoad)
					{
						Out[0] ++;
						Out[Out[0]] = m;
					}
				}
   			
				if (Add)
				{
					for (m = Out[0]; m > 0; m--)
					{
						for (l = Out[m]; l < ParetoSetSize; l++)
						{
							ParetoSetInsertions[l] = ParetoSetInsertions[l+1];
						}
						ParetoSetSize --;
					}
     			
					ParetoSetSize ++;
					ParetoSetInsertions[ParetoSetSize] = CandInsertions[NO_CandInsertions];
				}
  		    	  
				w = inst_tasks[LeftTasks[n]].inverse;

				NO_CandInsertions ++;
				CandInsertions[NO_CandInsertions].InsertedTask = w;
				CandInsertions[NO_CandInsertions].InsertRouteID = j;
				CandInsertions[NO_CandInsertions].InsertPos = k;
				CandInsertions[NO_CandInsertions].InsertCost = min_cost[inst_tasks[Routes1[j][k-1]].tail_node][inst_tasks[w].head_node]
				+min_cost[inst_tasks[w].tail_node][inst_tasks[Routes1[j][k]].head_node]
				-min_cost[inst_tasks[Routes1[j][k-1]].tail_node][inst_tasks[Routes1[j][k]].head_node];
		    
				CandInsertions[NO_CandInsertions].InsertVioLoad = IVLoad;
   			
				Out[0] = 0;
				Add = 1;
   			
				for (m = 1; m <= ParetoSetSize; m++)
				{
					if (CandInsertions[NO_CandInsertions].InsertCost > ParetoSetInsertions[m].InsertCost
						&& CandInsertions[NO_CandInsertions].InsertVioLoad > ParetoSetInsertions[m].InsertVioLoad)
					{
						Add = 0;
						break;
					}
					else if (CandInsertions[NO_CandInsertions].InsertCost < ParetoSetInsertions[m].InsertCost
						&& CandInsertions[NO_CandInsertions].InsertVioLoad < ParetoSetInsertions[m].InsertVioLoad)
					{
						Out[0] ++;
						Out[Out[0]] = m;
					}
				}
   			
				if (Add)
				{
					for (m = Out[0]; m > 0; m--)
					{
						for (l = Out[m]; l < ParetoSetSize; l++)
						{
							ParetoSetInsertions[l] = ParetoSetInsertions[l+1];
						}
						ParetoSetSize --;
					}
     			
					ParetoSetSize ++;
					ParetoSetInsertions[ParetoSetSize] = CandInsertions[NO_CandInsertions];
				}
			}
	  }
	  
	  //if (NO_CurrRoutes < vehicle_num)
	  //{
		  NO_CandInsertions ++;
		  CandInsertions[NO_CandInsertions].InsertedTask = LeftTasks[n];
		  CandInsertions[NO_CandInsertions].InsertRouteID = 0;
		  CandInsertions[NO_CandInsertions].InsertPos = 2;
		  CandInsertions[NO_CandInsertions].InsertCost = min_cost[DEPOT][inst_tasks[LeftTasks[n]].head_node]+
			  min_cost[inst_tasks[LeftTasks[n]].tail_node][DEPOT];
		  CandInsertions[NO_CandInsertions].InsertVioLoad = 0;
	  //}
  	  
	  Out[0] = 0;
	  Add = 1;

	  for (m = 1; m <= ParetoSetSize; m++)
	  {
		  if (CandInsertions[NO_CandInsertions].InsertCost > ParetoSetInsertions[m].InsertCost
			  && CandInsertions[NO_CandInsertions].InsertVioLoad > ParetoSetInsertions[m].InsertVioLoad)
		  {
			  Add = 0;
			  break;
		  }
		  else if (CandInsertions[NO_CandInsertions].InsertCost < ParetoSetInsertions[m].InsertCost
			  && CandInsertions[NO_CandInsertions].InsertVioLoad < ParetoSetInsertions[m].InsertVioLoad)
		  {
			  Out[0] ++;
			  Out[Out[0]] = m;
		  }
	  }
   		
	  if (Add)
	  {
		  for (m = Out[0]; m > 0; m--)
		  {
			  for (l = Out[m]; l < ParetoSetSize; l++)
			  {
				  ParetoSetInsertions[l] = ParetoSetInsertions[l+1];
			  }
			  ParetoSetSize --;
		  }
    			
		  ParetoSetSize ++;
		  ParetoSetInsertions[ParetoSetSize] = CandInsertions[NO_CandInsertions];
	  }
		
	  k = rand_choose(ParetoSetSize);
	  BestInsertion = ParetoSetInsertions[k];
		
	  if (BestInsertion.InsertRouteID == 0)
	  {
		  Routes1[0][0] ++;
		  Routes1[Routes1[0][0]][0] = 3;
		  Routes1[Routes1[0][0]][1] = 0;
		  Routes1[Routes1[0][0]][2] = BestInsertion.InsertedTask;
		  Routes1[Routes1[0][0]][3] = 0;
			
		  XCLds[0] ++;
		  XCLds[XCLds[0]] = inst_tasks[BestInsertion.InsertedTask].demand;
	  }
	  else
	  {
		  add_element(Routes1[BestInsertion.InsertRouteID], BestInsertion.InsertedTask, BestInsertion.InsertPos);
		  XCLds[BestInsertion.InsertRouteID] += inst_tasks[BestInsertion.InsertedTask].demand;
	  }
		
		/*for (i = 1; i <= LeftTasks[0]; i++)
		{
			if (LeftTasks[i] == BestInsertion.InsertedTask || LeftTasks[i] == inst_tasks[BestInsertion.InsertedTask].Inv)
				break;
		}
		delete_element(LeftTasks, i);*/
	}
	
	//if (Prt)
  //	printf("add ok\n");
	
	xed_child->sequence[0] = 1;
	for (i = 1; i <= Routes1[0][0]; i++)
	{
		if (Routes1[i][0] == 2)
			continue;
		
		xed_child->sequence[0] --;
		link_array(xed_child->sequence, Routes1[i]);
	}
	
	xed_child->total_cost = get_task_seq_total_cost(xed_child->sequence, inst_tasks);
	memcpy(xed_child->route_seg_load, XCLds, sizeof(XCLds));
	
	for (i = xed_child->route_seg_load[0]; i > 0; i--)
	{
		if (xed_child->route_seg_load[i] == 0)
			delete_element(xed_child->route_seg_load, i);
	}
	
	xed_child->total_vio_load = get_total_vio_load(xed_child->route_seg_load);
	//get_route_seg_length(xed_child->route_seg_length, xed_child->sequence, inst_tasks);
	//xed_child->max_length = max(xed_child->route_seg_length);
	
	int RouteLoad;
 	find_ele_positions(Positions, xed_child->sequence, 0);
 	for (i = 1; i < Positions[0]; i++)
 	{
 		RouteLoad = 0;
 		for (j = Positions[i]; j < Positions[i+1]; j++)
 		{
 			RouteLoad += inst_tasks[xed_child->sequence[j]].demand;
 		}
 		
 		if (RouteLoad != xed_child->route_seg_load[i])
 		{
 			printf("XChild Seq\n");
     	for (k = 1; k <= xed_child->sequence[0]; k++)
 	    {
     		printf("%d  ", xed_child->sequence[k]);
 	    }
     	printf("\n");
     	
     	printf("i = %d, RouteLoad = %d\n", i, RouteLoad);
     	
     	printf("XChild Loads\n");
     	for (k = 1; k <= xed_child->route_seg_load[0]; k++)
 	    {
     		printf("%d  ", xed_child->route_seg_load[k]);
 	    }
     	printf("\n");
     	exit(0);
 		}
 	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void lns_mut(individual *c, individual *p, const task *inst_tasks)
{
	indi_copy(c, p);
	double coef = best_fsb_solution.total_cost/capacity*(1.0*best_fsb_solution.total_cost/p->total_cost+1.0*p->total_vio_load/capacity+1);
	c->fitness = c->total_cost+coef*c->total_vio_load;

	int const_fsb = 0, const_infsb = 0;
	int count = 0;

	int imp = 1;
	while (imp)
	{
		count ++;
		imp = 0;

		if (c->total_vio_load == 0)
		{
			const_fsb ++;
		}
		else
		{
			const_infsb ++;
		}
		
		if (count%5 == 0)
		{
			if (const_fsb == 5)
			{
				coef /= 2;
				c->fitness = c->total_cost+coef*c->total_vio_load;
			}
			else if (const_infsb == 5)
			{
				coef *= 2;
				c->fitness = c->total_cost+coef*c->total_vio_load;
			}
			
			const_fsb = 0;
			const_infsb = 0;
		}

		individual tmp_indi;
		indi_copy(&tmp_indi, c);

		lns(c, coef, 1, inst_tasks);

		if (c->fitness < tmp_indi.fitness)
			imp = 1;

		if (c->total_vio_load == 0 && c->total_cost < best_fsb_solution.total_cost)
		{
			indi_copy(&best_fsb_solution, c);
			printf("new best cost = %d\n", best_fsb_solution.total_cost);
		}
	}

	imp = 1;
	while (imp)
	{
		count ++;
		imp = 0;

		if (c->total_vio_load == 0)
		{
			const_fsb ++;
		}
		else
		{
			const_infsb ++;
		}
		
		if (count%5 == 0)
		{
			if (const_fsb == 5)
			{
				coef /= 2;
				c->fitness = c->total_cost+coef*c->total_vio_load;
			}
			else if (const_infsb == 5)
			{
				coef *= 2;
				c->fitness = c->total_cost+coef*c->total_vio_load;
			}
			
			const_fsb = 0;
			const_infsb = 0;
		}

		individual tmp_indi;
		indi_copy(&tmp_indi, c);

		lns(c, coef, 2, inst_tasks);

		if (c->fitness < tmp_indi.fitness)
			imp = 1;

		if (c->total_vio_load == 0 && c->total_cost < best_fsb_solution.total_cost)
		{
			indi_copy(&best_fsb_solution, c);
			printf("new best cost = %d\n", best_fsb_solution.total_cost);
		}
	}

	imp = 1;
	while (imp)
	{
		count ++;
		imp = 0;

		if (c->total_vio_load == 0)
		{
			const_fsb ++;
		}
		else
		{
			const_infsb ++;
		}
		
		if (count%5 == 0)
		{
			if (const_fsb == 5)
			{
				coef /= 2;
				c->fitness = c->total_cost+coef*c->total_vio_load;
			}
			else if (const_infsb == 5)
			{
				coef *= 2;
				c->fitness = c->total_cost+coef*c->total_vio_load;
			}
			
			const_fsb = 0;
			const_infsb = 0;
		}

		individual tmp_indi;
		indi_copy(&tmp_indi, c);

		lns(c, coef, 1, inst_tasks);

		if (c->fitness < tmp_indi.fitness)
			imp = 1;

		if (c->total_vio_load == 0 && c->total_cost < best_fsb_solution.total_cost)
		{
			indi_copy(&best_fsb_solution, c);
			printf("new best cost = %d\n", best_fsb_solution.total_cost);
		}
	}

	/*if (c->total_vio_load == 0)
	{
		int positions[MAX_SEG_TAG_LENGTH];

		find_ele_positions(positions, c->sequence, 0);
		for (int i = 1; i < positions[0]; i++)
		{
			int route[MAX_TASK_ROUTE_LENGTH];
			int fh_route1[MAX_TASK_ROUTE_LENGTH], fh_route2[MAX_TASK_ROUTE_LENGTH], fh_route[MAX_TASK_ROUTE_LENGTH];
			copy_sub_array(route, c->sequence, positions[i], positions[i+1]);
			fred_heuristic(fh_route1, route, inst_tasks, 0);
			fred_heuristic(fh_route2, route, inst_tasks, 1);
			int cost1 = get_task_seq_total_cost(fh_route1, inst_tasks);
			int cost2 = get_task_seq_total_cost(fh_route2, inst_tasks);
			int cost3;
			if (cost1 < cost2)
			{
				memcpy(fh_route, fh_route1, sizeof(fh_route1));
				cost3 = cost1;
			}
			else
			{
				memcpy(fh_route, fh_route2, sizeof(fh_route2));
				cost3 = cost2;
			}
			int cost4 = get_task_seq_total_cost(route, inst_tasks);
			if (cost3 < cost4)
			{
				for (int j = positions[i]; j < positions[i+1]; j++)
				{
					c->sequence[j] = fh_route[j-positions[i]+1];
				}
				c->total_cost += cost3-cost4;
				c->route_seg_length[i] += cost3-cost4;
			}
		}
		c->max_length = max(c->route_seg_length);
	}*/
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void lns(individual *indi, double coef, int nsize, const task *inst_tasks)
{
	indi->fitness = indi->total_cost+coef*indi->total_vio_load;

	if (nsize == 1) // traditional move operators, i.e., single insertion, double insertion, swap, etc.
	{
		move si_move, di_move, swap_move, next_move;
		next_move.fitness = INF;

		single_insertion(&si_move, indi, coef, inst_tasks);
		double_insertion(&di_move, indi, coef, inst_tasks);
		swap(&swap_move, indi, coef, inst_tasks);

		if (si_move.fitness < next_move.fitness)
			next_move = si_move;
		if (di_move.fitness < next_move.fitness)
			next_move = di_move;
		if (swap_move.fitness < next_move.fitness)
			next_move = swap_move;

		//printf("next type = %d, task1 = %d, orig_seg = %d, tar_seg = %d, orig_pos = %d, tar_pos = %d, total_cost = %d, fitness = %lf\n", next_move.type, next_move.task1,
		//	next_move.orig_seg, next_move.targ_seg, next_move.orig_pos, next_move.targ_pos, next_move.total_cost, next_move.fitness);

		int orig_ptr, targ_ptr, seg_ptr1, seg_ptr2;
		orig_ptr = 0;
		targ_ptr = 0;
		seg_ptr1 = 0;
		seg_ptr2 = 0;
		for (int i = 1; i < indi->sequence[0]; i++)
		{
			if (indi->sequence[i] == 0)
			{
				if (seg_ptr1 < next_move.orig_seg)
					seg_ptr1 ++;
				if (seg_ptr2 < next_move.targ_seg)
					seg_ptr2 ++;
				if (seg_ptr1 == next_move.orig_seg && orig_ptr == 0)
					orig_ptr = i+next_move.orig_pos-1;
				if (seg_ptr2 == next_move.targ_seg && targ_ptr == 0)
					targ_ptr = i+next_move.targ_pos-1;
			}
			if (orig_ptr != 0 && targ_ptr != 0)
				break;
		}

		//printf("before\n");
		//print_one_dim_array(indi->sequence);
		//print_one_dim_array(indi->route_seg_length);
		//printf("totalcost = %d, maxlength = %d\n", indi->total_cost, indi->max_length);

		switch (next_move.type)
		{
		case SI:
			{
				delete_element(indi->sequence, orig_ptr);
				if (targ_ptr > orig_ptr)
					targ_ptr --;
				//indi->route_seg_length[next_move.orig_seg] = next_move.orig_length;
				indi->route_seg_load[next_move.orig_seg] -= inst_tasks[next_move.task1].demand;
				if (next_move.targ_seg > indi->route_seg_load[0])
				{
					indi->sequence[0] ++;
					indi->sequence[indi->sequence[0]] = next_move.task1;
					indi->sequence[0] ++;
					indi->sequence[indi->sequence[0]] = 0;
					//indi->route_seg_length[0] ++;
					//indi->route_seg_length[indi->route_seg_length[0]] = next_move.targ_length;
					indi->route_seg_load[0] ++;
					indi->route_seg_load[indi->route_seg_load[0]] = inst_tasks[next_move.task1].demand;
				}
				else
				{
					add_element(indi->sequence, next_move.task1, targ_ptr);
					//indi->route_seg_length[next_move.targ_seg] = next_move.targ_length;
					indi->route_seg_load[next_move.targ_seg] += inst_tasks[next_move.task1].demand;
				}
			}
			break;
		case DI:
			{
				delete_element(indi->sequence, orig_ptr+1);
				delete_element(indi->sequence, orig_ptr);
				if (targ_ptr > orig_ptr)
					targ_ptr -= 2;
				//indi->route_seg_length[next_move.orig_seg] = next_move.orig_length;
				indi->route_seg_load[next_move.orig_seg] -= inst_tasks[next_move.task1].demand+inst_tasks[next_move.task2].demand;
				if (next_move.targ_seg > indi->route_seg_load[0])
				{
					indi->sequence[0] ++;
					indi->sequence[indi->sequence[0]] = next_move.task1;
					indi->sequence[0] ++;
					indi->sequence[indi->sequence[0]] = next_move.task2;
					indi->sequence[0] ++;
					indi->sequence[indi->sequence[0]] = 0;
					//indi->route_seg_length[0] ++;
					//indi->route_seg_length[indi->route_seg_length[0]] = next_move.targ_length;
					indi->route_seg_load[0] ++;
					indi->route_seg_load[indi->route_seg_load[0]] = inst_tasks[next_move.task1].demand+inst_tasks[next_move.task2].demand;
				}
				else
				{
					add_element(indi->sequence, next_move.task2, targ_ptr);
					add_element(indi->sequence, next_move.task1, targ_ptr);
					//indi->route_seg_length[next_move.targ_seg] = next_move.targ_length;
					indi->route_seg_load[next_move.targ_seg] += inst_tasks[next_move.task1].demand+inst_tasks[next_move.task2].demand;
				}
			}
			break;
		case SWAP:
			{
				indi->sequence[targ_ptr] = next_move.task1;
				indi->sequence[orig_ptr] = next_move.task2;
				//indi->route_seg_length[next_move.orig_seg] = next_move.orig_length;
				//indi->route_seg_length[next_move.targ_seg] = next_move.targ_length;
				indi->route_seg_load[next_move.orig_seg] -= inst_tasks[next_move.task1].demand-inst_tasks[next_move.task2].demand;
				indi->route_seg_load[next_move.targ_seg] += inst_tasks[next_move.task1].demand-inst_tasks[next_move.task2].demand;
			}
			break;
		}
		//indi->max_length = next_move.max_length;
		indi->total_cost = next_move.total_cost;
		indi->total_vio_load = next_move.total_vio_load;
		indi->fitness = next_move.fitness;

		if (indi->route_seg_load[next_move.orig_seg] == 0)
		{
			if (next_move.type == DI && next_move.orig_seg > next_move.targ_seg)
			{
				delete_element(indi->sequence, orig_ptr+1);
			}
			else
			{
				delete_element(indi->sequence, orig_ptr);
			}
			//delete_element(indi->route_seg_length, next_move.orig_seg);
			delete_element(indi->route_seg_load, next_move.orig_seg);
		}

		//printf("after\n");
		//print_one_dim_array(indi->sequence);
		//print_one_dim_array(indi->route_seg_length);
		//printf("totalcost = %d, maxlength = %d\n", indi->total_cost, indi->max_length);
	}
	else
	{
		individual tmp_indi, next_indi;
		next_indi.fitness = INF;
		//int task_routes[MAX_SEG_TAG_LENGTH][MAX_TASK_SEG_LENGTH];

		task_routes[0][0] = 1;
		task_routes[1][0] = 1;
		task_routes[1][1] = 0;
		for (int i = 2; i <= indi->sequence[0]; i++)
		{
			task_routes[task_routes[0][0]][0] ++;
			task_routes[task_routes[0][0]][task_routes[task_routes[0][0]][0]] = indi->sequence[i];

			if (indi->sequence[i] == 0 && i < indi->sequence[0])
			{
				task_routes[0][0] ++;
				task_routes[task_routes[0][0]][0] = 1;
				task_routes[task_routes[0][0]][1] = 0;
			}
		}

		if (task_routes[0][0] < nsize)
			return;
		
		int multi = task_routes[0][0];
		long long int ub_trial = task_routes[0][0];
		for (int i = 1; i < nsize; i++)
		{
			multi --;
			ub_trial *= multi;
		}		
		multi = nsize;
		for (int i = 1; i < nsize; i++)
		{
			ub_trial /= multi;
			multi --;
		}
		
		int maxcount = ub_trial;
		if (maxcount > MAX_ENSSIZE)
			maxcount = MAX_ENSSIZE;

		typedef struct lns_comb
		{
			int ids[MAX_NSIZE+1];
		}lns_comb;
		
		lns_comb cand_combs[MAX_ENSSIZE+1];
		
		int pointers[MAX_NSIZE+1];
		for (int i = 1; i <= nsize; i++)
		{
			pointers[i] = i;
		}
		
		int curr_ptr;
		for (int i = 1; i <= maxcount; i++)
		{
			cand_combs[i].ids[0] = nsize;
			for (int j = 1; j <= nsize; j++)
			{
				cand_combs[i].ids[j] = pointers[j];
			}
			
			curr_ptr = nsize;
			while (pointers[curr_ptr] == task_routes[0][0]-nsize+curr_ptr)
			{
				curr_ptr --;
			}
			
			if (curr_ptr == 0)
				break;
			
			pointers[curr_ptr] ++;
			for (int j = curr_ptr+1; j <= nsize; j++)
			{
				pointers[j] = pointers[j-1]+1;
			}
		}

		int lns_routes[MAX_NSIZE+1];
		for (int i = 0; i < maxcount; i++)
		{
			memcpy(lns_routes, cand_combs[i].ids, sizeof(cand_combs[i].ids));

			int sel_total_load = 0;
			for (int j = 1; j <= lns_routes[0]; j++)
			{
				sel_total_load += indi->route_seg_load[lns_routes[j]];
			}

			if (sel_total_load > nsize*capacity)
				continue;

			int serve_mark[MAX_TASK_TAG_LENGTH];
			memset(serve_mark, 0, sizeof(serve_mark));
			serve_mark[0] = task_num;
			for (int j = 1; j <= lns_routes[0]; j++)
			{
				for (int k = 2; k < task_routes[lns_routes[j]][0]; k++)
				{
					serve_mark[task_routes[lns_routes[j]][k]] = 1;
					serve_mark[inst_tasks[task_routes[lns_routes[j]][k]].inverse] = 1;
				}
			}

			path_scanning(&tmp_indi, inst_tasks, serve_mark);

			for (int j = 1; j <= task_routes[0][0]; j++)
			{
				int lnsed = 0;
				for (int k = 1; k <= lns_routes[0]; k++)
				{
					if (j == lns_routes[k])
					{
						lnsed = 1;
						break;
					}
				}
				if (lnsed)
					continue;

				tmp_indi.sequence[0] --;
				link_array(tmp_indi.sequence, task_routes[j]);
				tmp_indi.route_seg_load[0] ++;
				tmp_indi.route_seg_load[tmp_indi.route_seg_load[0]] = indi->route_seg_load[j];
				if (indi->route_seg_load[j] > capacity)
					tmp_indi.total_vio_load += indi->route_seg_load[j];
			}

			tmp_indi.total_cost = get_task_seq_total_cost(tmp_indi.sequence, inst_tasks);
			//get_route_seg_length(tmp_indi.route_seg_length, tmp_indi.sequence, inst_tasks);
			//tmp_indi.max_length = max(tmp_indi.route_seg_length);
			tmp_indi.fitness = tmp_indi.total_cost+coef*tmp_indi.total_vio_load;

			if (tmp_indi.fitness < next_indi.fitness)
				indi_copy(&next_indi, &tmp_indi);
		}

		/*long long int ub_trial = 2;
		for (int i = 1; i < task_routes[0][0]; i++)
		{
			ub_trial *= 2;
		}

		int lns_routes[MAX_SEG_TAG_LENGTH];
		int maxcount = 100;
		int count = 0;
		for (int i = 0; i < ub_trial; i++)
		{
			lns_routes[0] = 0;
			for (int j = 1; j <= task_routes[0][0]; j++)
			{
				if (bit_one(i, j))
				{
					lns_routes[0] ++;
					lns_routes[lns_routes[0]] = j;
				}
			}

			if (lns_routes[0] != nsize)
				continue;

			count ++;

			int sel_total_load = 0;
			for (int j = 1; j <= lns_routes[0]; j++)
			{
				sel_total_load += indi->route_seg_load[lns_routes[j]];
			}

			if (sel_total_load > nsize*capacity)
				continue;

			int serve_mark[MAX_TASK_TAG_LENGTH];
			memset(serve_mark, 0, sizeof(serve_mark));
			serve_mark[0] = task_num;
			for (int j = 1; j <= lns_routes[0]; j++)
			{
				for (int k = 2; k < task_routes[lns_routes[j]][0]; k++)
				{
					serve_mark[task_routes[lns_routes[j]][k]] = 1;
					serve_mark[inst_tasks[task_routes[lns_routes[j]][k]].inverse] = 1;
				}
			}

			path_scanning(&tmp_indi, inst_tasks, serve_mark);

			for (int j = 1; j <= task_routes[0][0]; j++)
			{
				if (bit_one(i, j))
					continue;

				tmp_indi.sequence[0] --;
				link_array(tmp_indi.sequence, task_routes[j]);
				tmp_indi.route_seg_load[0] ++;
				tmp_indi.route_seg_load[tmp_indi.route_seg_load[0]] = indi->route_seg_load[j];
				if (indi->route_seg_load[j] > capacity)
					tmp_indi.total_vio_load += indi->route_seg_load[j];
			}

			tmp_indi.total_cost = get_task_seq_total_cost(tmp_indi.sequence, inst_tasks);
			get_route_seg_length(tmp_indi.route_seg_length, tmp_indi.sequence, inst_tasks);
			tmp_indi.max_length = max(tmp_indi.route_seg_length);
			tmp_indi.fitness = tmp_indi.total_cost+coef*tmp_indi.total_vio_load;

			if (tmp_indi.fitness < next_indi.fitness)
				indi_copy(&next_indi, &tmp_indi);

			if (count == maxcount)
				break;
		}*/

		if (next_indi.fitness < indi->fitness)
			indi_copy(indi, &next_indi);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void single_insertion(move *best_move, individual *indi, double coef, const task *inst_tasks)
{
	best_move->type = SI;
	best_move->fitness = INF;

	//int tmp_route_seg_length[MAX_SEG_TAG_LENGTH];

	
	task_routes[0][0] = 1;
	task_routes[1][0] = 1;
	task_routes[1][1] = 0;
	for (int i = 2; i <= indi->sequence[0]; i++)
	{
		task_routes[task_routes[0][0]][0] ++;
		task_routes[task_routes[0][0]][task_routes[task_routes[0][0]][0]] = indi->sequence[i];

		if (indi->sequence[i] == 0 && i < indi->sequence[0])
		{
			task_routes[0][0] ++;
			task_routes[task_routes[0][0]][0] = 1;
			task_routes[task_routes[0][0]][1] = 0;
		}
	}

	move tmp_move;
	for (int s1 = 1; s1 <= task_routes[0][0]; s1++)
	{
		tmp_move.orig_seg = s1;
		for (int i = 2; i < task_routes[s1][0]; i++)
		{
			tmp_move.orig_pos = i;
			for (int s2 = 1; s2 <= task_routes[0][0]+1; s2++) /* s2 > task_routes[0][0] --> create a new route */
			{
				if (s2 == s1)
					continue;

				tmp_move.targ_seg = s2;
				if (s2 > task_routes[0][0]/* && task_routes[0][0] < vehicle_num*/)
				{
					tmp_move.targ_pos = 2;
					tmp_move.total_vio_load = indi->total_vio_load;
					if (indi->route_seg_load[s1] > capacity)
						tmp_move.total_vio_load -= indi->route_seg_load[s1]-capacity;
					if (indi->route_seg_load[s1]-inst_tasks[task_routes[s1][i]].demand > capacity)
						tmp_move.total_vio_load += indi->route_seg_load[s1]-inst_tasks[task_routes[s1][i]].demand-capacity;

					tmp_move.task1 = task_routes[s1][i];
					//memcpy(tmp_route_seg_length, indi->route_seg_length, sizeof(indi->route_seg_length));
					//tmp_route_seg_length[s1] +=
					//	min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
					//	min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
					//	inst_tasks[task_routes[s1][i]].serv_cost-
					//	min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node];
					//tmp_route_seg_length[0] ++;
					//tmp_route_seg_length[tmp_route_seg_length[0]] =
					//	min_cost[DEPOT][inst_tasks[tmp_move.task1].head_node]+inst_tasks[tmp_move.task1].serv_cost+
					//	min_cost[inst_tasks[tmp_move.task1].tail_node][DEPOT];

					//tmp_move.max_length = max(tmp_route_seg_length);

					//tmp_move.orig_length = tmp_route_seg_length[s1];
					//tmp_move.targ_length = tmp_route_seg_length[tmp_route_seg_length[0]];
					//tmp_move.total_cost = indi->total_cost+
					//	tmp_route_seg_length[s1]+tmp_route_seg_length[tmp_route_seg_length[0]]-indi->route_seg_length[s1];

					tmp_move.total_cost = indi->total_cost+
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]+
						min_cost[DEPOT][inst_tasks[tmp_move.task1].head_node]+
						min_cost[inst_tasks[tmp_move.task1].tail_node][DEPOT]-
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node];

					tmp_move.fitness = tmp_move.total_cost+coef*tmp_move.total_vio_load;

					if (tmp_move.fitness < best_move->fitness)
					{
						best_move->task1 = tmp_move.task1;
						best_move->orig_seg = tmp_move.orig_seg;
						best_move->targ_seg = tmp_move.targ_seg;
						best_move->orig_pos = tmp_move.orig_pos;
						best_move->targ_pos = tmp_move.targ_pos;
						//best_move->orig_length = tmp_move.orig_length;
						//best_move->targ_length = tmp_move.targ_length;
						best_move->total_cost = tmp_move.total_cost;
						best_move->total_vio_load = tmp_move.total_vio_load;
						//best_move->max_length = tmp_move.max_length;
						best_move->fitness = tmp_move.fitness;
					}
					continue;
				}

				for (int j = 2; j <= task_routes[s2][0]; j++)
				{
					if (inst_tasks[task_routes[s2][j-1]].tail_node == inst_tasks[task_routes[s2][j]].head_node)
						continue;

					tmp_move.targ_pos = j;
					tmp_move.total_vio_load = indi->total_vio_load;
					if (indi->route_seg_load[s1] > capacity)
						tmp_move.total_vio_load -= indi->route_seg_load[s1]-capacity;
					if (indi->route_seg_load[s2] > capacity)
						tmp_move.total_vio_load -= indi->route_seg_load[s2]-capacity;
					if (indi->route_seg_load[s1]-inst_tasks[task_routes[s1][i]].demand > capacity)
						tmp_move.total_vio_load += indi->route_seg_load[s1]-inst_tasks[task_routes[s1][i]].demand-capacity;
					if (indi->route_seg_load[s2]+inst_tasks[task_routes[s1][i]].demand > capacity)
						tmp_move.total_vio_load += indi->route_seg_load[s2]+inst_tasks[task_routes[s1][i]].demand-capacity;

					tmp_move.task1 = task_routes[s1][i];
					//memcpy(tmp_route_seg_length, indi->route_seg_length, sizeof(indi->route_seg_length));
					//tmp_route_seg_length[s1] += 
					//	min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
					//	min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
					//	inst_tasks[task_routes[s1][i]].serv_cost-
					//	min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node];
					//tmp_route_seg_length[s2] +=
					//	min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
					//	inst_tasks[tmp_move.task1].serv_cost+
					//	min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
					//	min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node];
						
					//tmp_move.max_length = max(tmp_route_seg_length);

					//tmp_move.orig_length = tmp_route_seg_length[s1];
					//tmp_move.targ_length = tmp_route_seg_length[s2];
					//tmp_move.total_cost = indi->total_cost+
					//	tmp_route_seg_length[s1]+tmp_route_seg_length[s2]-indi->route_seg_length[s1]-indi->route_seg_length[s2];

					tmp_move.total_cost = indi->total_cost+
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]+
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
						min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node];

					tmp_move.fitness = tmp_move.total_cost+coef*tmp_move.total_vio_load;

					if (tmp_move.fitness < best_move->fitness)
					{
						best_move->task1 = tmp_move.task1;
						best_move->orig_seg = tmp_move.orig_seg;
						best_move->targ_seg = tmp_move.targ_seg;
						best_move->orig_pos = tmp_move.orig_pos;
						best_move->targ_pos = tmp_move.targ_pos;
						//best_move->orig_length = tmp_move.orig_length;
						//best_move->targ_length = tmp_move.targ_length;
						best_move->total_cost = tmp_move.total_cost;
						best_move->total_vio_load = tmp_move.total_vio_load;
						//best_move->max_length = tmp_move.max_length;
						best_move->fitness = tmp_move.fitness;
					}

					tmp_move.task1 = inst_tasks[task_routes[s1][i]].inverse;
					//memcpy(tmp_route_seg_length, indi->route_seg_length, sizeof(indi->route_seg_length));
					//tmp_route_seg_length[s1] += 
					//	min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
					//	min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
					//	inst_tasks[task_routes[s1][i]].serv_cost-
					//	min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node];
					//tmp_route_seg_length[s2] +=
					//	min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
					//	inst_tasks[tmp_move.task1].serv_cost+
					//	min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
					//	min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node];
						
					//tmp_move.max_length = max(tmp_route_seg_length);

					//tmp_move.orig_length = tmp_route_seg_length[s1];
					//tmp_move.targ_length = tmp_route_seg_length[s2];
					//tmp_move.total_cost = indi->total_cost+
					//	tmp_route_seg_length[s1]+tmp_route_seg_length[s2]-indi->route_seg_length[s1]-indi->route_seg_length[s2];

					tmp_move.total_cost = indi->total_cost+
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]+
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
						min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node];

					tmp_move.fitness = tmp_move.total_cost+coef*tmp_move.total_vio_load;

					if (tmp_move.fitness < best_move->fitness)
					{
						best_move->task1 = tmp_move.task1;
						best_move->orig_seg = tmp_move.orig_seg;
						best_move->targ_seg = tmp_move.targ_seg;
						best_move->orig_pos = tmp_move.orig_pos;
						best_move->targ_pos = tmp_move.targ_pos;
						//best_move->orig_length = tmp_move.orig_length;
						//best_move->targ_length = tmp_move.targ_length;
						best_move->total_cost = tmp_move.total_cost;
						best_move->total_vio_load = tmp_move.total_vio_load;
						//best_move->max_length = tmp_move.max_length;
						best_move->fitness = tmp_move.fitness;
					}
				}
			}
		}
	}

	//printf("task = %d, orig_route = %d, targ_route = %d, orig_load = %d, tar_load = %d, add_cost = %d, add_vioload = %d, add_fit = %lf\n", best_move->task1,
	//	best_move->orig_route, best_move->targ_route, best_move->orig_load, best_move->targ_load, best_move->add_cost, best_move->total_vio_load,
	//	best_move->add_fit);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void double_insertion(move *best_move, individual *indi, double coef, const task *inst_tasks)
{
	best_move->type = DI;
	best_move->fitness = INF;

	//int tmp_route_seg_length[MAX_SEG_TAG_LENGTH];

	//int task_routes[MAX_SEG_TAG_LENGTH][MAX_TASK_SEG_LENGTH];
	task_routes[0][0] = 1;
	task_routes[1][0] = 1;
	task_routes[1][1] = 0;
	for (int i = 2; i <= indi->sequence[0]; i++)
	{
		task_routes[task_routes[0][0]][0] ++;
		task_routes[task_routes[0][0]][task_routes[task_routes[0][0]][0]] = indi->sequence[i];

		if (indi->sequence[i] == 0 && i < indi->sequence[0])
		{
			task_routes[0][0] ++;
			task_routes[task_routes[0][0]][0] = 1;
			task_routes[task_routes[0][0]][1] = 0;
		}
	}

	move tmp_move;
	for (int s1 = 1; s1 <= task_routes[0][0]; s1++)
	{
		if (task_routes[s1][0] < 4)
			continue;

		tmp_move.orig_seg = s1;
		for (int i = 2; i < task_routes[s1][0]-1; i++)
		{
			tmp_move.orig_pos = i;
			for (int s2 = 1; s2 <= task_routes[0][0]+1; s2++) /* s2 > task_routes[0][0] --> create a new route */
			{
				if (s2 == s1)
					continue;

				tmp_move.targ_seg = s2;
				if (s2 > task_routes[0][0]/* && task_routes[0][0] < vehicle_num*/)
				{
					if (task_routes[s1][0] <= 4)
						continue;

					tmp_move.targ_pos = 2;
					tmp_move.total_vio_load = indi->total_vio_load;
					if (indi->route_seg_load[s1] > capacity)
						tmp_move.total_vio_load -= indi->route_seg_load[s1]-capacity;
					if (indi->route_seg_load[s1]-inst_tasks[task_routes[s1][i]].demand-inst_tasks[task_routes[s1][i+1]].demand > capacity)
						tmp_move.total_vio_load += indi->route_seg_load[s1]-inst_tasks[task_routes[s1][i]].demand-inst_tasks[task_routes[s1][i+1]].demand-capacity;

					tmp_move.task1 = task_routes[s1][i];
					tmp_move.task2 = task_routes[s1][i+1];
					//memcpy(tmp_route_seg_length, indi->route_seg_length, sizeof(indi->route_seg_length));
					//tmp_route_seg_length[s1] +=
					//	min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]-
					//	min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
					//	inst_tasks[task_routes[s1][i]].serv_cost-
					//	min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
					//	inst_tasks[task_routes[s1][i+1]].serv_cost-
					//	min_cost[inst_tasks[task_routes[s1][i+1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node];
					//tmp_route_seg_length[0] ++;
					//tmp_route_seg_length[tmp_route_seg_length[0]] =
					//	min_cost[DEPOT][inst_tasks[tmp_move.task1].head_node]+inst_tasks[tmp_move.task1].serv_cost+
					//	min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[tmp_move.task2].head_node]+inst_tasks[tmp_move.task2].serv_cost+
					//	min_cost[inst_tasks[tmp_move.task2].tail_node][DEPOT];

					//tmp_move.max_length = max(tmp_route_seg_length);

					//tmp_move.orig_length = tmp_route_seg_length[s1];
					//tmp_move.targ_length = tmp_route_seg_length[tmp_route_seg_length[0]];
					//tmp_move.total_cost = indi->total_cost+
					//	tmp_route_seg_length[s1]+tmp_route_seg_length[tmp_route_seg_length[0]]-indi->route_seg_length[s1];

					tmp_move.total_cost = indi->total_cost+
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i+1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]+
						min_cost[DEPOT][inst_tasks[tmp_move.task1].head_node]+
						min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[tmp_move.task2].head_node]+
						min_cost[inst_tasks[tmp_move.task2].tail_node][DEPOT];

					tmp_move.fitness = tmp_move.total_cost+coef*tmp_move.total_vio_load;

					if (tmp_move.fitness < best_move->fitness)
					{
						best_move->task1 = tmp_move.task1;
						best_move->task2 = tmp_move.task2;
						best_move->orig_seg = tmp_move.orig_seg;
						best_move->targ_seg = tmp_move.targ_seg;
						best_move->orig_pos = tmp_move.orig_pos;
						best_move->targ_pos = tmp_move.targ_pos;
						//best_move->orig_length = tmp_move.orig_length;
						//best_move->targ_length = tmp_move.targ_length;
						best_move->total_cost = tmp_move.total_cost;
						best_move->total_vio_load = tmp_move.total_vio_load;
						//best_move->max_length = tmp_move.max_length;
						best_move->fitness = tmp_move.fitness;
					}

					tmp_move.task1 = task_routes[s1][i];
					tmp_move.task2 = inst_tasks[task_routes[s1][i+1]].inverse;
					//memcpy(tmp_route_seg_length, indi->route_seg_length, sizeof(indi->route_seg_length));
					//tmp_route_seg_length[s1] +=
					//	min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]-
					//	min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
					//	inst_tasks[task_routes[s1][i]].serv_cost-
					//	min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
					//	inst_tasks[task_routes[s1][i+1]].serv_cost-
					//	min_cost[inst_tasks[task_routes[s1][i+1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node];
					//tmp_route_seg_length[0] ++;
					//tmp_route_seg_length[tmp_route_seg_length[0]] =
					//	min_cost[DEPOT][inst_tasks[tmp_move.task1].head_node]+inst_tasks[tmp_move.task1].serv_cost+
					//	min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[tmp_move.task2].head_node]+inst_tasks[tmp_move.task2].serv_cost+
					//	min_cost[inst_tasks[tmp_move.task2].tail_node][DEPOT];

					//tmp_move.max_length = max(tmp_route_seg_length);

					//tmp_move.orig_length = tmp_route_seg_length[s1];
					//tmp_move.targ_length = tmp_route_seg_length[tmp_route_seg_length[0]];
					//tmp_move.total_cost = indi->total_cost+
					//	tmp_route_seg_length[s1]+tmp_route_seg_length[tmp_route_seg_length[0]]-indi->route_seg_length[s1];

					tmp_move.total_cost = indi->total_cost+
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i+1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]+
						min_cost[DEPOT][inst_tasks[tmp_move.task1].head_node]+
						min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[tmp_move.task2].head_node]+
						min_cost[inst_tasks[tmp_move.task2].tail_node][DEPOT];

					tmp_move.fitness = tmp_move.total_cost+coef*tmp_move.total_vio_load;

					if (tmp_move.fitness < best_move->fitness)
					{
						best_move->task1 = tmp_move.task1;
						best_move->task2 = tmp_move.task2;
						best_move->orig_seg = tmp_move.orig_seg;
						best_move->targ_seg = tmp_move.targ_seg;
						best_move->orig_pos = tmp_move.orig_pos;
						best_move->targ_pos = tmp_move.targ_pos;
						//best_move->orig_length = tmp_move.orig_length;
						//best_move->targ_length = tmp_move.targ_length;
						best_move->total_cost = tmp_move.total_cost;
						best_move->total_vio_load = tmp_move.total_vio_load;
						//best_move->max_length = tmp_move.max_length;
						best_move->fitness = tmp_move.fitness;
					}
					continue;
				}

				for (int j = 2; j <= task_routes[s2][0]; j++) 
				{
					if (inst_tasks[task_routes[s2][j-1]].tail_node == inst_tasks[task_routes[s2][j]].head_node)
						continue;

					tmp_move.targ_pos = j;
					tmp_move.total_vio_load = indi->total_vio_load;
					if (indi->route_seg_load[s1] > capacity)
						tmp_move.total_vio_load -= indi->route_seg_load[s1]-capacity;
					if (indi->route_seg_load[s2] > capacity)
						tmp_move.total_vio_load -= indi->route_seg_load[s2]-capacity;
					if (indi->route_seg_load[s1]-inst_tasks[task_routes[s1][i]].demand-inst_tasks[task_routes[s1][i+1]].demand > capacity)
						tmp_move.total_vio_load += indi->route_seg_load[s1]-inst_tasks[task_routes[s1][i]].demand-inst_tasks[task_routes[s1][i+1]].demand-capacity;
					if (indi->route_seg_load[s2]+inst_tasks[task_routes[s1][i]].demand+inst_tasks[task_routes[s1][i+1]].demand > capacity)
						tmp_move.total_vio_load += indi->route_seg_load[s2]+inst_tasks[task_routes[s1][i]].demand+inst_tasks[task_routes[s1][i+1]].demand-capacity;

					tmp_move.task1 = task_routes[s1][i];
					tmp_move.task2 = task_routes[s1][i+1];
					//memcpy(tmp_route_seg_length, indi->route_seg_length, sizeof(indi->route_seg_length));
					//tmp_route_seg_length[s1] += 
					//	min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]-
					//	min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
					//	inst_tasks[task_routes[s1][i]].serv_cost-
					//	min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
					//	inst_tasks[task_routes[s1][i+1]].serv_cost-
					//	min_cost[inst_tasks[task_routes[s1][i+1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node];
					//tmp_route_seg_length[s2] +=
					//	min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
					//	inst_tasks[tmp_move.task1].serv_cost+
					//	min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[tmp_move.task2].head_node]+
					//	inst_tasks[tmp_move.task2].serv_cost+
					//	min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
					//	min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node];
						
					//tmp_move.max_length = max(tmp_route_seg_length);

					//tmp_move.orig_length = tmp_route_seg_length[s1];
					//tmp_move.targ_length = tmp_route_seg_length[s2];
					//tmp_move.total_cost = indi->total_cost+
					//	tmp_route_seg_length[s1]+tmp_route_seg_length[s2]-indi->route_seg_length[s1]-indi->route_seg_length[s2];

					tmp_move.total_cost = indi->total_cost+
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i+1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]+
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
						min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[tmp_move.task2].head_node]+
						min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node];

					tmp_move.fitness = tmp_move.total_cost+coef*tmp_move.total_vio_load;

					if (tmp_move.fitness < best_move->fitness)
					{
						best_move->task1 = tmp_move.task1;
						best_move->task2 = tmp_move.task2;
						best_move->orig_seg = tmp_move.orig_seg;
						best_move->targ_seg = tmp_move.targ_seg;
						best_move->orig_pos = tmp_move.orig_pos;
						best_move->targ_pos = tmp_move.targ_pos;
						//best_move->orig_length = tmp_move.orig_length;
						//best_move->targ_length = tmp_move.targ_length;
						best_move->total_cost = tmp_move.total_cost;
						best_move->total_vio_load = tmp_move.total_vio_load;
						//best_move->max_length = tmp_move.max_length;
						best_move->fitness = tmp_move.fitness;
					}

					////

					tmp_move.task1 = inst_tasks[task_routes[s1][i]].inverse;
					tmp_move.task2 = task_routes[s1][i+1];
					//memcpy(tmp_route_seg_length, indi->route_seg_length, sizeof(indi->route_seg_length));
					//tmp_route_seg_length[s1] += 
					//	min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]-
					//	min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
					//	inst_tasks[task_routes[s1][i]].serv_cost-
					//	min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
					//	inst_tasks[task_routes[s1][i+1]].serv_cost-
					//	min_cost[inst_tasks[task_routes[s1][i+1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node];
					//tmp_route_seg_length[s2] +=
					//	min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
					//	inst_tasks[tmp_move.task1].serv_cost+
					//	min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[tmp_move.task2].head_node]+
					//	inst_tasks[tmp_move.task2].serv_cost+
					//	min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
					//	min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node];
						
					//tmp_move.max_length = max(tmp_route_seg_length);

					//tmp_move.orig_length = tmp_route_seg_length[s1];
					//tmp_move.targ_length = tmp_route_seg_length[s2];
					//tmp_move.total_cost = indi->total_cost+
					//	tmp_route_seg_length[s1]+tmp_route_seg_length[s2]-indi->route_seg_length[s1]-indi->route_seg_length[s2];

					tmp_move.total_cost = indi->total_cost+
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i+1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]+
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
						min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[tmp_move.task2].head_node]+
						min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node];

					tmp_move.fitness = tmp_move.total_cost+coef*tmp_move.total_vio_load;

					if (tmp_move.fitness < best_move->fitness)
					{
						best_move->task1 = tmp_move.task1;
						best_move->task2 = tmp_move.task2;
						best_move->orig_seg = tmp_move.orig_seg;
						best_move->targ_seg = tmp_move.targ_seg;
						best_move->orig_pos = tmp_move.orig_pos;
						best_move->targ_pos = tmp_move.targ_pos;
						//best_move->orig_length = tmp_move.orig_length;
						//best_move->targ_length = tmp_move.targ_length;
						best_move->total_cost = tmp_move.total_cost;
						best_move->total_vio_load = tmp_move.total_vio_load;
						//best_move->max_length = tmp_move.max_length;
						best_move->fitness = tmp_move.fitness;
					}

					////

					tmp_move.task1 = task_routes[s1][i];
					tmp_move.task2 = inst_tasks[task_routes[s1][i+1]].inverse;
					//memcpy(tmp_route_seg_length, indi->route_seg_length, sizeof(indi->route_seg_length));
					//tmp_route_seg_length[s1] += 
					//	min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]-
					//	min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
					//	inst_tasks[task_routes[s1][i]].serv_cost-
					//	min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
					//	inst_tasks[task_routes[s1][i+1]].serv_cost-
					//	min_cost[inst_tasks[task_routes[s1][i+1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node];
					//tmp_route_seg_length[s2] +=
					//	min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
					//	inst_tasks[tmp_move.task1].serv_cost+
					//	min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[tmp_move.task2].head_node]+
					//	inst_tasks[tmp_move.task2].serv_cost+
					//	min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
					//	min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node];
						
					//tmp_move.max_length = max(tmp_route_seg_length);

					//tmp_move.orig_length = tmp_route_seg_length[s1];
					//tmp_move.targ_length = tmp_route_seg_length[s2];
					//tmp_move.total_cost = indi->total_cost+
					//	tmp_route_seg_length[s1]+tmp_route_seg_length[s2]-indi->route_seg_length[s1]-indi->route_seg_length[s2];

					tmp_move.total_cost = indi->total_cost+
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i+1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]+
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
						min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[tmp_move.task2].head_node]+
						min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node];

					tmp_move.fitness = tmp_move.total_cost+coef*tmp_move.total_vio_load;

					if (tmp_move.fitness < best_move->fitness)
					{
						best_move->task1 = tmp_move.task1;
						best_move->task2 = tmp_move.task2;
						best_move->orig_seg = tmp_move.orig_seg;
						best_move->targ_seg = tmp_move.targ_seg;
						best_move->orig_pos = tmp_move.orig_pos;
						best_move->targ_pos = tmp_move.targ_pos;
						//best_move->orig_length = tmp_move.orig_length;
						//best_move->targ_length = tmp_move.targ_length;
						best_move->total_cost = tmp_move.total_cost;
						best_move->total_vio_load = tmp_move.total_vio_load;
						//best_move->max_length = tmp_move.max_length;
						best_move->fitness = tmp_move.fitness;
					}

					////

					tmp_move.task1 = inst_tasks[task_routes[s1][i]].inverse;
					tmp_move.task2 = inst_tasks[task_routes[s1][i+1]].inverse;
					//memcpy(tmp_route_seg_length, indi->route_seg_length, sizeof(indi->route_seg_length));
					//tmp_route_seg_length[s1] += 
					//	min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]-
					//	min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
					//	inst_tasks[task_routes[s1][i]].serv_cost-
					//	min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
					//	inst_tasks[task_routes[s1][i+1]].serv_cost-
					//	min_cost[inst_tasks[task_routes[s1][i+1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node];
					//tmp_route_seg_length[s2] +=
					//	min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
					//	inst_tasks[tmp_move.task1].serv_cost+
					//	min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[tmp_move.task2].head_node]+
					//	inst_tasks[tmp_move.task2].serv_cost+
					//	min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
					//	min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node];
						
					//tmp_move.max_length = max(tmp_route_seg_length);

					//tmp_move.orig_length = tmp_route_seg_length[s1];
					//tmp_move.targ_length = tmp_route_seg_length[s2];
					//tmp_move.total_cost = indi->total_cost+
					//	tmp_route_seg_length[s1]+tmp_route_seg_length[s2]-indi->route_seg_length[s1]-indi->route_seg_length[s2];

					tmp_move.total_cost = indi->total_cost+
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i+1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]+
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
						min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[tmp_move.task2].head_node]+
						min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node];

					tmp_move.fitness = tmp_move.total_cost+coef*tmp_move.total_vio_load;

					if (tmp_move.fitness < best_move->fitness)
					{
						best_move->task1 = tmp_move.task1;
						best_move->task2 = tmp_move.task2;
						best_move->orig_seg = tmp_move.orig_seg;
						best_move->targ_seg = tmp_move.targ_seg;
						best_move->orig_pos = tmp_move.orig_pos;
						best_move->targ_pos = tmp_move.targ_pos;
						//best_move->orig_length = tmp_move.orig_length;
						//best_move->targ_length = tmp_move.targ_length;
						best_move->total_cost = tmp_move.total_cost;
						best_move->total_vio_load = tmp_move.total_vio_load;
						//best_move->max_length = tmp_move.max_length;
						best_move->fitness = tmp_move.fitness;
					}
				}
			}
		}
	}

	//printf("task1 = %d, task2 = %d, orig_route = %d, targ_route = %d, orig_load = %d, tar_load = %d, add_cost = %d, add_vioload = %d, add_fit = %lf\n",
	//	best_move->task1, best_move->task2, best_move->orig_route, best_move->targ_route, best_move->orig_load, best_move->targ_load, best_move->add_cost,
	//	best_move->total_vio_load, best_move->add_fit);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void swap(move *best_move, individual *indi, double coef, const task *inst_tasks)
{
	best_move->type = SWAP;
	best_move->fitness = INF;

	//int tmp_route_seg_length[MAX_SEG_TAG_LENGTH];

	//int task_routes[MAX_SEG_TAG_LENGTH][MAX_TASK_SEG_LENGTH];
	task_routes[0][0] = 1;
	task_routes[1][0] = 1;
	task_routes[1][1] = 0;
	for (int i = 2; i <= indi->sequence[0]; i++)
	{
		task_routes[task_routes[0][0]][0] ++;
		task_routes[task_routes[0][0]][task_routes[task_routes[0][0]][0]] = indi->sequence[i];

		if (indi->sequence[i] == 0 && i < indi->sequence[0])
		{
			task_routes[0][0] ++;
			task_routes[task_routes[0][0]][0] = 1;
			task_routes[task_routes[0][0]][1] = 0;
		}
	}

	move tmp_move;
	for (int s1 = 1; s1 < task_routes[0][0]; s1++)
	{
		tmp_move.orig_seg = s1;
		for (int i = 2; i < task_routes[s1][0]; i++)
		{
			if (inst_tasks[task_routes[s1][i-1]].tail_node == inst_tasks[task_routes[s1][i]].head_node &&
				inst_tasks[task_routes[s1][i]].tail_node == inst_tasks[task_routes[s1][i+1]].head_node)
				continue;

			tmp_move.orig_pos = i;
			for (int s2 = s1+1; s2 <= task_routes[0][0]; s2++)
			{
				tmp_move.targ_seg = s2;
				for (int j = 2; j < task_routes[s2][0]; j++)
				{
					if (inst_tasks[task_routes[s2][j-1]].tail_node == inst_tasks[task_routes[s2][j]].head_node &&
						inst_tasks[task_routes[s2][j]].tail_node == inst_tasks[task_routes[s2][j+1]].head_node)
						continue;

					tmp_move.targ_pos = j;

					tmp_move.total_vio_load = indi->total_vio_load;
					if (indi->route_seg_load[s1] > capacity)
						tmp_move.total_vio_load -= indi->route_seg_load[s1]-capacity;
					if (indi->route_seg_load[s2] > capacity)
						tmp_move.total_vio_load -= indi->route_seg_load[s2]-capacity;
					if (indi->route_seg_load[s1]-inst_tasks[task_routes[s1][i]].demand+inst_tasks[task_routes[s2][j]].demand > capacity)
						tmp_move.total_vio_load += indi->route_seg_load[s1]-inst_tasks[task_routes[s1][i]].demand+inst_tasks[task_routes[s2][j]].demand-capacity;
					if (indi->route_seg_load[s2]+inst_tasks[task_routes[s1][i]].demand-inst_tasks[task_routes[s2][j]].demand > capacity)
						tmp_move.total_vio_load += indi->route_seg_load[s2]+inst_tasks[task_routes[s1][i]].demand-inst_tasks[task_routes[s2][j]].demand-capacity;

					tmp_move.task1 = task_routes[s1][i];
					tmp_move.task2 = task_routes[s2][j];
					//memcpy(tmp_route_seg_length, indi->route_seg_length, sizeof(indi->route_seg_length));
					//tmp_route_seg_length[s1] += 
					//	min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[tmp_move.task2].head_node]+
					//	inst_tasks[tmp_move.task2].serv_cost+
					//	min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
					//	min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
					//	inst_tasks[task_routes[s1][i]].serv_cost-
					//	min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node];
					//tmp_route_seg_length[s2] +=
					//	min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
					//	inst_tasks[tmp_move.task1].serv_cost+
					//	min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[task_routes[s2][j+1]].head_node]-
					//	min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
					//	inst_tasks[task_routes[s2][j]].serv_cost-
					//	min_cost[inst_tasks[task_routes[s2][j]].tail_node][inst_tasks[task_routes[s2][j+1]].head_node];
						
					//tmp_move.max_length = max(tmp_route_seg_length);

					//tmp_move.orig_length = tmp_route_seg_length[s1];
					//tmp_move.targ_length = tmp_route_seg_length[s2];
					//tmp_move.total_cost = indi->total_cost+
					//	tmp_route_seg_length[s1]+tmp_route_seg_length[s2]-indi->route_seg_length[s1]-indi->route_seg_length[s2];

					tmp_move.total_cost = indi->total_cost+
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[tmp_move.task2].head_node]+
						min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]+
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
						min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[task_routes[s2][j+1]].head_node]-
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
						min_cost[inst_tasks[task_routes[s2][j]].tail_node][inst_tasks[task_routes[s2][j+1]].head_node];

					tmp_move.fitness = tmp_move.total_cost+coef*tmp_move.total_vio_load;

					if (tmp_move.fitness < best_move->fitness)
					{
						best_move->task1 = tmp_move.task1;
						best_move->task2 = tmp_move.task2;
						best_move->orig_seg = tmp_move.orig_seg;
						best_move->targ_seg = tmp_move.targ_seg;
						best_move->orig_pos = tmp_move.orig_pos;
						best_move->targ_pos = tmp_move.targ_pos;
						//best_move->orig_length = tmp_move.orig_length;
						//best_move->targ_length = tmp_move.targ_length;
						best_move->total_cost = tmp_move.total_cost;
						best_move->total_vio_load = tmp_move.total_vio_load;
						//best_move->max_length = tmp_move.max_length;
						best_move->fitness = tmp_move.fitness;
					}

					////

					tmp_move.task1 = inst_tasks[task_routes[s1][i]].inverse;
					tmp_move.task2 = task_routes[s2][j];
					//memcpy(tmp_route_seg_length, indi->route_seg_length, sizeof(indi->route_seg_length));
					//tmp_route_seg_length[s1] += 
					//	min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[tmp_move.task2].head_node]+
					//	inst_tasks[tmp_move.task2].serv_cost+
					//	min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
					//	min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
					//	inst_tasks[task_routes[s1][i]].serv_cost-
					//	min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node];
					//tmp_route_seg_length[s2] +=
					//	min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
					//	inst_tasks[tmp_move.task1].serv_cost+
					//	min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[task_routes[s2][j+1]].head_node]-
					//	min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
					//	inst_tasks[task_routes[s2][j]].serv_cost-
					//	min_cost[inst_tasks[task_routes[s2][j]].tail_node][inst_tasks[task_routes[s2][j+1]].head_node];
						
					//tmp_move.max_length = max(tmp_route_seg_length);

					//tmp_move.orig_length = tmp_route_seg_length[s1];
					//tmp_move.targ_length = tmp_route_seg_length[s2];
					//tmp_move.total_cost = indi->total_cost+
					//	tmp_route_seg_length[s1]+tmp_route_seg_length[s2]-indi->route_seg_length[s1]-indi->route_seg_length[s2];

					tmp_move.total_cost = indi->total_cost+
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[tmp_move.task2].head_node]+
						min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]+
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
						min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[task_routes[s2][j+1]].head_node]-
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
						min_cost[inst_tasks[task_routes[s2][j]].tail_node][inst_tasks[task_routes[s2][j+1]].head_node];

					tmp_move.fitness = tmp_move.total_cost+coef*tmp_move.total_vio_load;

					if (tmp_move.fitness < best_move->fitness)
					{
						best_move->task1 = tmp_move.task1;
						best_move->task2 = tmp_move.task2;
						best_move->orig_seg = tmp_move.orig_seg;
						best_move->targ_seg = tmp_move.targ_seg;
						best_move->orig_pos = tmp_move.orig_pos;
						best_move->targ_pos = tmp_move.targ_pos;
						//best_move->orig_length = tmp_move.orig_length;
						//best_move->targ_length = tmp_move.targ_length;
						best_move->total_cost = tmp_move.total_cost;
						best_move->total_vio_load = tmp_move.total_vio_load;
						//best_move->max_length = tmp_move.max_length;
						best_move->fitness = tmp_move.fitness;
					}

					////

					tmp_move.task1 = task_routes[s1][i];
					tmp_move.task2 = inst_tasks[task_routes[s2][j]].inverse;
					//memcpy(tmp_route_seg_length, indi->route_seg_length, sizeof(indi->route_seg_length));
					//tmp_route_seg_length[s1] += 
					//	min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[tmp_move.task2].head_node]+
					//	inst_tasks[tmp_move.task2].serv_cost+
					//	min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
					//	min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
					//	inst_tasks[task_routes[s1][i]].serv_cost-
					//	min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node];
					//tmp_route_seg_length[s2] +=
					//	min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
					//	inst_tasks[tmp_move.task1].serv_cost+
					//	min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[task_routes[s2][j+1]].head_node]-
					//	min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
					//	inst_tasks[task_routes[s2][j]].serv_cost-
					//	min_cost[inst_tasks[task_routes[s2][j]].tail_node][inst_tasks[task_routes[s2][j+1]].head_node];
						
					//tmp_move.max_length = max(tmp_route_seg_length);

					//tmp_move.orig_length = tmp_route_seg_length[s1];
					//tmp_move.targ_length = tmp_route_seg_length[s2];
					//tmp_move.total_cost = indi->total_cost+
					//	tmp_route_seg_length[s1]+tmp_route_seg_length[s2]-indi->route_seg_length[s1]-indi->route_seg_length[s2];

					tmp_move.total_cost = indi->total_cost+
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[tmp_move.task2].head_node]+
						min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]+
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
						min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[task_routes[s2][j+1]].head_node]-
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
						min_cost[inst_tasks[task_routes[s2][j]].tail_node][inst_tasks[task_routes[s2][j+1]].head_node];

					tmp_move.fitness = tmp_move.total_cost+coef*tmp_move.total_vio_load;

					if (tmp_move.fitness < best_move->fitness)
					{
						best_move->task1 = tmp_move.task1;
						best_move->task2 = tmp_move.task2;
						best_move->orig_seg = tmp_move.orig_seg;
						best_move->targ_seg = tmp_move.targ_seg;
						best_move->orig_pos = tmp_move.orig_pos;
						best_move->targ_pos = tmp_move.targ_pos;
						//best_move->orig_length = tmp_move.orig_length;
						//best_move->targ_length = tmp_move.targ_length;
						best_move->total_cost = tmp_move.total_cost;
						best_move->total_vio_load = tmp_move.total_vio_load;
						//best_move->max_length = tmp_move.max_length;
						best_move->fitness = tmp_move.fitness;
					}

					////

					tmp_move.task1 = inst_tasks[task_routes[s1][i]].inverse;
					tmp_move.task2 = inst_tasks[task_routes[s2][j]].inverse;
					//memcpy(tmp_route_seg_length, indi->route_seg_length, sizeof(indi->route_seg_length));
					//tmp_route_seg_length[s1] += 
					//	min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[tmp_move.task2].head_node]+
					//	inst_tasks[tmp_move.task2].serv_cost+
					//	min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
					//	min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
					//	inst_tasks[task_routes[s1][i]].serv_cost-
					//	min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node];
					//tmp_route_seg_length[s2] +=
					//	min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
					//	inst_tasks[tmp_move.task1].serv_cost+
					//	min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[task_routes[s2][j+1]].head_node]-
					//	min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
					//	inst_tasks[task_routes[s2][j]].serv_cost-
					//	min_cost[inst_tasks[task_routes[s2][j]].tail_node][inst_tasks[task_routes[s2][j+1]].head_node];
						
					//tmp_move.max_length = max(tmp_route_seg_length);

					//tmp_move.orig_length = tmp_route_seg_length[s1];
					//tmp_move.targ_length = tmp_route_seg_length[s2];
					//tmp_move.total_cost = indi->total_cost+
					//	tmp_route_seg_length[s1]+tmp_route_seg_length[s2]-indi->route_seg_length[s1]-indi->route_seg_length[s2];

					tmp_move.total_cost = indi->total_cost+
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[tmp_move.task2].head_node]+
						min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]+
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
						min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[task_routes[s2][j+1]].head_node]-
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
						min_cost[inst_tasks[task_routes[s2][j]].tail_node][inst_tasks[task_routes[s2][j+1]].head_node];

					tmp_move.fitness = tmp_move.total_cost+coef*tmp_move.total_vio_load;

					if (tmp_move.fitness < best_move->fitness)
					{
						best_move->task1 = tmp_move.task1;
						best_move->task2 = tmp_move.task2;
						best_move->orig_seg = tmp_move.orig_seg;
						best_move->targ_seg = tmp_move.targ_seg;
						best_move->orig_pos = tmp_move.orig_pos;
						best_move->targ_pos = tmp_move.targ_pos;
						//best_move->orig_length = tmp_move.orig_length;
						//best_move->targ_length = tmp_move.targ_length;
						best_move->total_cost = tmp_move.total_cost;
						best_move->total_vio_load = tmp_move.total_vio_load;
						//best_move->max_length = tmp_move.max_length;
						best_move->fitness = tmp_move.fitness;
					}
				}
			}
		}
	}

	//printf("task1 = %d, task2 = %d, orig_route = %d, targ_route = %d, orig_load = %d, tar_load = %d, add_cost = %d, add_vioload = %d, add_fit = %lf\n",
	//	best_move->task1, best_move->task2, best_move->orig_route, best_move->targ_route, best_move->orig_load, best_move->targ_load, best_move->add_cost,
	//	best_move->total_vio_load, best_move->add_fit);
}