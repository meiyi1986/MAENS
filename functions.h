
#define INF 2100000000
#define DUMMY_CYCLE 0
#define MAX_EDGE_TASK_NUM 500
#define MAX_ARC_TASK_NUM 0
#define MAX_TASK_NUM 1000
#define MAX_NODE_NUM 300
#define MAX_ARCS_TAG_LENGTH 1001
#define MAX_TASK_TAG_LENGTH 1001
#define MAX_NODE_TAG_LENGTH 300
#define MAX_ROUTE_TAG_LENGTH 50
#define MAX_SEG_TAG_LENGTH 50
#define MAX_TASK_SEG_LENGTH 550
#define MAX_NODE_SEG_LENGTH 1000
#define MAX_TASK_ROUTE_LENGTH 550
#define MAX_NODE_ROUTE_LENGTH 1000
#define MAX_TASK_SEQ_LENGTH 550
#define MAX_FACS_TAG_LENGTH 5
#define SI 1
#define DI 2
#define SWAP 3

#define MAX_INST_NUM 50

#define MAX_POPSIZE 30
#define MAX_TOTALSIZE 210
#define MAX_NONDOMINATED_NUM 1000

#define M_trial 10
#define M_PROB 0.2
#define M_ite 100
#define M_wite 100

#define MAX_NSIZE 10 // upper bound of nsize
#define MAX_ENSSIZE 100 // maximal ENS neighborhood size

extern int vertex_num;
extern int req_edge_num;
extern int req_arc_num;
extern int nonreq_edge_num;
extern int nonreq_arc_num;
extern int task_num;
extern int total_arc_num;
extern int vehicle_num;
extern int capacity;
extern int lower_bound;

extern int DEPOT;

extern int LBs[MAX_INST_NUM];
extern int inst_LB;

extern int trav_cost[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
extern int serve_cost[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
extern int shortest_path[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
extern int min_cost[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];

/* randomness */

int* rand_perm(int num);
// generate a random permutation from 1 to num

int rand_choose(int num);
// choose a number randomly between 1 and num (num must be not too large)

/* operations for arrays */

void print_one_dim_array(int *a);

void print_two_dim_matrix(int **a, int row, int col);

void delete_element(int *a, int k);
// delete the kth element from array a

void add_element(int *a, int e, int k);
// add element e to the kth position of array a

void reverse_direction(int *a, int k1, int k2);
// reverse the direction of array a from position k1 to k2

void find_ele_positions(int *positions, int *a, int e);
// find the positions that element e occurs in array a

void copy_sub_array(int *b, int *a, int k1, int k2);
// copy sub-array a[k1:k2] to array b

void insert_array(int *b, int *a, int k);
// insert array a into position k of array b

void replace_sub_array(int *b, int *a, int k1, int k2);
// replace the sub-array b[k1:k2] with array a

void link_array(int *b, int *a);
// link array a behind array b

int equal(int *a, int *b);
// if array a equals array b, return 1, otherwise return 0

int min(int *a);
// get minimal value within array a

int max(int *a);
// get maximal value within array a

int set_bit(int num, int k);
// set kth bit of num to 1

int unset_bit(int num, int k);
// set kth bit of num to 0

int bit_one(int num, int k);
// return 1 if the kth bit of num is 1, otherwise return 0

/* operations for task sequences (a special array) */

typedef struct task
{
	int head_node;
	int tail_node;
	int dead_cost;
	int serv_cost;
	int demand;
	int inverse;
} task;

typedef struct arc
{
	int tail_node;
	int head_node;
	int trav_cost;
} arc;

typedef struct individual
{
	int sequence[MAX_TASK_SEQ_LENGTH];
	int route_seg_load[MAX_SEG_TAG_LENGTH];
	int total_cost;
	int total_vio_load;
	double fitness;
	int pred[MAX_TASK_TAG_LENGTH];
	int succ[MAX_TASK_TAG_LENGTH];
} individual;

extern individual pop[MAX_TOTALSIZE];
extern individual best_fsb_solution;
extern int popsize;

int delete_element(individual *indis, int indis_size, int k);
// delete the kth individual from the individual set indis with its size as indis_size

int find_arc(int head_node, int tail_node, const arc *inst_arcs);
// find arc index according to head_node and tail_node

int find_task(int head_node, int tail_node, const task *inst_tasks);
// find task index according to head_node and tail_node

void remove_task_seq_delimiters(int *task_seq);
// remove delimiters '0' from task_seq

int split(int *split_task_seq, int *one_task_seq, const task *inst_tasks);
// split operator

int fleet_limited_split(int *split_task_seq, int *one_task_seq, const task *inst_tasks);
// split operator subject to vehicle_num

int legal(int *task_seq, const task *inst_tasks);
// if task_seq is legal, return 1, otherwise return 0

int get_task_seq_total_cost(int *task_seq, const task *inst_tasks);
// get total cost of task_seq

void get_route_loads(int *route_loads, int *task_seq, const task *inst_tasks);
// get load of each route according to task_seq

int get_total_vio_load(int *route_loads);
// get total_vio_load according to route_loads

void get_route_seg_length(int *route_seg_length, int *task_seq, const task *inst_tasks);
// get length of each route according to task_seq

void get_assignment(int *assignment, int *task_seq, const task *inst_tasks);
// get assignment of each task according to task_seq

void get_predsucc(int *pred, int *succ, int *task_seq, const task *inst_tasks);
// get predecessors and successors of the task sequence

int calc_distance(individual *indi1, individual *indi2, const task *inst_tasks);
// calculate the distance between indi1 and indi2 according to their 'pred' and 'succ'

void indi_copy(individual *target, individual *source);
// copy source to target

/* initialization and preprocessing functions */

void mod_dijkstra();

void path_scanning(individual *ps_indi, const task *inst_tasks, int *serve_mark);

void rand_scanning(individual *rs_indi, const task *inst_tasks, int *serve_mark);

/* search operators */

typedef struct move
{
	int type;
	int task1;
	int task2;
	int orig_seg;
	int targ_seg;
	int orig_pos;
	int targ_pos;
	int orig_length;
	int targ_length;
	int total_cost;
	int total_vio_load;
	int max_length;
	double fitness;
} move;

void rand_selection(int &id1, int &id2, individual *pop);

void tour_selection(int &id1, int &id2, individual *pop);

void SBX(individual *xed_child, individual *p1, individual *p2, const task *inst_tasks);

void lns_mut(individual *c, individual *p, const task *inst_tasks);

void lns(individual *indi, double coef, int nsize, const task *inst_tasks);

void single_insertion(move *best_move, individual *curr_indi, double coef, const task *inst_tasks);

void double_insertion(move *best_move, individual *curr_indi, double coef, const task *inst_tasks);

void swap(move *best_move, individual *curr_indi, double coef, const task *inst_tasks);

void global_repair_operator(individual *indi, const task *inst_tasks);