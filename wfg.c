// To do: 
// - can we sort less often or reduce/optimise dominance checks? 
// - should we use FPL's data structure? 
// - two changes in read.c 
// - heuristics 

// opt:  0 = basic, 1 = sorting, 2 = slicing to 2D, 3 = slicing to 3D 
// mode: 0 = hv,    1 = exclhv,  2 = smallest point complete, 3 = smallest point inorder, 4 = smallest point BFQ

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "wfg.h"
#include "avl.h"

int n;     // the number of objectives 
POINT ref; // the reference point 

#define BEATS(x,y)   (x <  y) // change these for max/minimisation
#define BEATSEQ(x,y) (x <= y) // currently minimisation 

#define WORSE(x,y)   (BEATS(y,x) ? (x) : (y)) 
#define BETTER(x,y)  (BEATS(y,x) ? (y) : (x)) 

#define BIG          BETTER(999999, -999999)

#define mode2        (mode == 0 ? 1 : 0) // used in memory management 
#define REFERENCE_POINT 100.0

FRONT *fs;      // memory management stuff 
int fr = 0;     // current depth 
int frmax = -1; // max depth malloced so far (for opt = 0) 
int maxm = 0;   // identify the biggest fronts in the file 
int maxn = 0;


static avl_tree_t *tree;
double hv(FRONT);

static int compare_tree_asc( const void *p1, const void *p2)
{
    const double x1= *((const double *)p1+1);
    const double x2= *((const double *)p2+1);

    if (x1 != x2)
        return (x1 > x2) ? -1 : 1;
    else
        return 0;
}



int greater(const void *v1, const void *v2)
// this sorts points improving in the last objective
{
  POINT p = *(POINT*)v1;
  POINT q = *(POINT*)v2;
  #if mode == 0 && opt == 1
  for (int i = n - fr - 1; i >= 0; i--)
  #else
  for (int i = n - 1; i >= 0; i--)
  #endif
    if BEATS(p.objectives[i],q.objectives[i]) return  1;
    else
    // can we manage with only one comparison?
    if BEATS(q.objectives[i],p.objectives[i]) return -1;
  return 0;
}


int smallersofar(const void *v1, const void *v2)
// this sorts jobs by partial volume 
{
  JOB j1 = *(JOB*)v1;
  JOB j2 = *(JOB*)v2;
  return j1.partial > j2.partial;
}


int dominates2way(POINT p, POINT q)
// returns -1 if p dominates q, 1 if q dominates p, 2 if p == q, 0 otherwise 
{
  int z = 2;
  // domination could be checked in either order 
  #if mode == 0 && opt == 1
  for (int i = n - fr - 1; i >= 0; i--)
  #else
  for (int i = n - 1; i >= 0; i--)
  #endif
    if BEATS(p.objectives[i],q.objectives[i]) {if (z ==  1) return 0; else z = -1;}
    else
    // can we manage with only one comparison?
    if BEATS(q.objectives[i],p.objectives[i]) {if (z == -1) return 0; else z =  1;}
  return z;
}


void makeDominatedBit(FRONT ps, int p)
// creates the front ps[p+1 ..] in fs[fr], with each point bounded by ps[p] and dominated points removed 
{
  // when opt = 0 each new frame is allocated as needed, because the worst-case needs #frames = #points 
  #if opt == 0
  if (fr > frmax)
    {frmax = fr;
     fs[fr].points = malloc(sizeof(POINT) * maxm);
     for (int j = 0; j < maxm; j++) 
     {
       fs[fr].points[j].objectives = malloc(sizeof(OBJECTIVE) * maxn);
     }
    }
  #endif

  int z = ps.nPoints - 1 - p;
  for (int i = 0; i < z; i++)
    for (int j = 0; j < n; j++) 
      fs[fr].points[i].objectives[j] = WORSE(ps.points[p].objectives[j],ps.points[p + 1 + i].objectives[j]); 
  POINT t;
  fs[fr].nPoints = 1;
  for (int i = 1; i < z; i++)
    {int j = 0;
     bool keep = true;
     while (j < fs[fr].nPoints && keep)
       switch (dominates2way(fs[fr].points[i], fs[fr].points[j]))
	 {case -1: t = fs[fr].points[j];
                   fs[fr].points[j] = fs[fr].points[fs[fr].nPoints - 1]; 
                   fs[fr].points[fs[fr].nPoints - 1] = t; 
                   fs[fr].nPoints--; 
                   break;
          case  0: j++; break;
          // case  2: printf("Identical points!\n");
	  default: keep = false;
	 }
     if (keep) {t = fs[fr].points[fs[fr].nPoints]; 
                fs[fr].points[fs[fr].nPoints] = fs[fr].points[i]; 
                fs[fr].points[i] = t; 
                fs[fr].nPoints++;}
    }
  fr++;
}


double hv2(FRONT ps)
// returns the hypervolume of ps[0 ..] in 2D 
// assumes that ps is sorted improving
{
  double volume = fabs((ps.points[0].objectives[0] - ref.objectives[0]) * 
                       (ps.points[0].objectives[1] - ref.objectives[1])); 
  for (int i = 1; i < ps.nPoints; i++) 
    volume += fabs((ps.points[i].objectives[0] - ref.objectives[0]) * 
                   (ps.points[i].objectives[1] - ps.points[i - 1].objectives[1]));
  return volume;
}



double hv3_AVL(FRONT ps)
// returns the hypervolume of ps[0 ..] in 3D 
// assumes that ps is sorted improving
{
  avl_init_node(ps.points[ps.nPoints-1].tnode,ps.points[ps.nPoints-1].objectives);
  avl_insert_top(tree,ps.points[ps.nPoints-1].tnode);

  double hypera = (ref.objectives[0] - ps.points[ps.nPoints-1].objectives[0]) *
    (ref.objectives[1] - ps.points[ps.nPoints-1].objectives[1]);

  double height;
  if (ps.nPoints == 1)
    height = ref.objectives[2] - ps.points[ps.nPoints-1].objectives[2];
  else
    height = ps.points[ps.nPoints-2].objectives[2] - ps.points[ps.nPoints-1].objectives[2];

  double hyperv = hypera * height;

  for (int i = ps.nPoints - 2; i >= 0; i--)
  {
    if (i == 0)
      height = ref.objectives[2] - ps.points[i].objectives[2];
    else
      height = ps.points[i-1].objectives[2] - ps.points[i].objectives[2];

      // search tree for point q to the right of current point
      const double * prv_ip, * nxt_ip;
      avl_node_t *tnode;

      avl_init_node(ps.points[i].tnode, ps.points[i].objectives);

      if (avl_search_closest(tree, ps.points[i].objectives, &tnode) <= 0) {
          nxt_ip = (double *)(tnode->item);
          tnode = tnode->prev;
      } else {
          nxt_ip = (tnode->next!=NULL)
              ? (double *)(tnode->next->item)
              : ref.objectives;
      }
                // if p is not dominated
                if (nxt_ip[0] > ps.points[i].objectives[0]) {

                  // insert p in tree
                    avl_insert_after(tree, tnode, ps.points[i].tnode);

                    if (tnode !=NULL) {
                        prv_ip = (double *)(tnode->item);

                        if (prv_ip[0] > ps.points[i].objectives[0]) {
                            const double * cur_ip;

                            tnode = ps.points[i].tnode->prev;
                            // cur_ip = point dominated by pp with highest [0]-coordinate
                            cur_ip = (double *)(tnode->item);

                            // for each point in s in tree dominated by p
                            while (tnode->prev) {
                                prv_ip = (double *)(tnode->prev->item);
                                // decrease area by contribution of s
                                hypera -= (prv_ip[1] - cur_ip[1])*(nxt_ip[0] - cur_ip[0]);
                                if (prv_ip[0] < ps.points[i].objectives[0])
                                    break; // prv is not dominated by pp
                                cur_ip = prv_ip;
                                // remove s from tree
                                avl_unlink_node(tree,tnode);
                                tnode = tnode->prev;
                            }

                            // remove s from tree
                            avl_unlink_node(tree,tnode);

                            if (!tnode->prev) {
                                // decrease area by contribution of s
                                hypera -= (ref.objectives[1] - cur_ip[1])*(nxt_ip[0] - cur_ip[0]);
                                prv_ip = ref.objectives;
                            }
                        }
                    } else
                        prv_ip = ref.objectives;

                    // increase area by contribution of p
                    hypera += (prv_ip[1] -
                        ps.points[i].objectives[1])*(nxt_ip[0] -
                          ps.points[i].objectives[0]);

                }

                if (height > 0)
                    hyperv += hypera * height;
        }
        avl_clear_tree(tree);

        return hyperv;

  }



double hv3(FRONT ps)
// returns the hypervolume of ps[0 ..] in 3D 
// assumes that ps is sorted improving
{
  POINT *base = fs[fr].points;
  base[0].objectives[0] = ref.objectives[0];
  base[0].objectives[1] = BIG; // BETTER(ps.points[ps.nPoints - 1].objectives[1] + 1,ps.points[ps.nPoints - 1].objectives[1] - 1); 
  base[1].objectives[0] = ps.points[ps.nPoints - 1].objectives[0];
  base[1].objectives[1] = ps.points[ps.nPoints - 1].objectives[1];
  base[2].objectives[0] = BIG; // BETTER(ps.points[0].objectives[0] + 1,ps.points[0].objectives[0] - 1); 
  base[2].objectives[1] = ref.objectives[1];
  int z = 2;
  double area = fabs((base[1].objectives[0] - base[0].objectives[0]) * 
                     (base[1].objectives[1] - base[2].objectives[1]));
  double volume = 0;
  for (int i = ps.nPoints - 2; i >= 0; i--)
    {volume += area * fabs(ps.points[i].objectives[2] - ps.points[i + 1].objectives[2]);
     int jl = 1; int jr = z; 
     while (jr > jl + 1)
       {int m = (jl + jr) / 2;
        if BEATS(ps.points[i].objectives[0],base[m].objectives[0]) jl = m; else jr = m;
       }
     if BEATS(base[jl].objectives[0],ps.points[i].objectives[0]) jr = jl; else jl = jr;
     while BEATSEQ(ps.points[i].objectives[1],base[jl - 1].objectives[1])
       {area -= fabs((base[jl - 1].objectives[0] - base[jl - 2].objectives[0]) * 
                     (base[jl - 1].objectives[1] - base[jr].objectives[1]));
        jl--;
       }
     switch (jr - jl)
       {case  0: z++;
                 for (int k = z; k > jl; k--) 
                   {base[k].objectives[0] = base[k - 1].objectives[0];
                    base[k].objectives[1] = base[k - 1].objectives[1];
                   }
                 break;
        case  1: break;
        default: z -= jr - jl - 1;
                 for (int k = jl; k < z; k++) 
                   {base[k + 1].objectives[0] = base[k + jr - jl].objectives[0];
                    base[k + 1].objectives[1] = base[k + jr - jl].objectives[1];
	           } 
       }
     base[jl].objectives[0] = ps.points[i].objectives[0];
     base[jl].objectives[1] = ps.points[i].objectives[1];
     area += fabs((base[jl].objectives[0] - base[jl - 1].objectives[0]) * 
                  (base[jl].objectives[1] - base[jl + 1].objectives[1]));
    }
  volume += area * fabs(ref.objectives[2] - ps.points[0].objectives[2]);
  return volume;
}




double inclhv(POINT p)
// returns the inclusive hypervolume of p
{
  double volume = 1;
  for (int i = 0; i < n; i++) 
    volume *= fabs(p.objectives[i] - ref.objectives[i]);
  return volume;
}


double exclhv(FRONT ps, int p)
// returns the exclusive hypervolume of ps[p] relative to ps[p+1 ..] 
{
  double volume = inclhv(ps.points[p]);
  if (ps.nPoints > p + 1) 
    {
     makeDominatedBit(ps, p);
     volume -= hv(fs[fr - 1]);
     fr--;
    }
  return volume;
}


int smallestpoint0(FRONT ps)
// returns the least contributor in ps[0 ..] 
// evaluates each point completely 
{
  int smallest = 0;
  double smallestvolume = exclhv(ps, 0);
  POINT t;
  for (int i = 1; i < ps.nPoints; i++)
    {t = ps.points[0];
     ps.points[0] = ps.points[i];
     ps.points[i] = t;
     double volume = exclhv(ps, 0);
     printf ("Volume %d, %.7e\n", i, volume);
     if (volume < smallestvolume) {smallest = i; smallestvolume = volume;}
    }
  return smallest;
}



void fixHeap(JOB *js, int next)
// fixes the heap js, next being the new smallest element
{
  /*
        JOB jz = js[0];
        js[0] = js[next];
        js[next] = jz;
        if (js[0].left == ps.nPoints)
          {js[0].left = next; js[next].left = ps.nPoints; js[next].right = ps.nPoints;}
        else if (js[0].right == ps.nPoints)
          {js[0].right = next; js[next].left = ps.nPoints; js[next].right = ps.nPoints;}
        else ;
*/
}


int smallestpoint1(FRONT ps)
// returns the least contributor in ps[0 ..] 
// evaluates the first point completely, then as little as possible of each subsequent point 
{
  int smallest = 0;
  double smallestvolume = exclhv(ps, 0);
  POINT t;
  for (int i = 1; i < ps.nPoints; i++)
    {t = ps.points[0];
     ps.points[0] = ps.points[i];
     ps.points[i] = t;
     makeDominatedBit(ps, 0);
     double volume = inclhv(ps.points[0]);
     #if opt > 0
     qsort(fs[0].points, fs[0].nPoints, sizeof(POINT), greater);
     #endif
     for (int j = 0; j < fs[0].nPoints; j++)
       volume -= inclhv(fs[0].points[j]);
     int j = fs[0].nPoints - 2; 
     while (volume < smallestvolume && j >= 0)
       {makeDominatedBit(fs[0], j);
        volume += hv(fs[1]);
        fr--;
        j--;
       }
     fr--;
     // printf("j(%d) = %d\n", i, j);
     if (volume < smallestvolume) {smallest = i; smallestvolume = volume;}
    }
  return smallest;
}



int smallestpoint2(FRONT ps)
// returns the least contributor in ps[0 ..] 
// evaluates each point until it's positive, then (repeatedly) as little as possible of the current smallest point 
{
  // this space should be allocated once, up to maxm
  JOB *js = malloc(sizeof(JOB) * (ps.nPoints + 1));
  POINT t;
  for (int i = 0; i < ps.nPoints; i++)
    {t = ps.points[0];
     ps.points[0] = ps.points[i];
     ps.points[i] = t;
     makeDominatedBit(ps, 0);
     // printf("size(%d) = %d\n", i, fs[0].nPoints);
     // generalise makeDominatedBit so the target front is a parameter, and save the copying?
     js[i].sprime.nPoints = fs[0].nPoints;
     js[i].sprime.points = malloc(sizeof(POINT) * js[i].sprime.nPoints);
     for (int j = 0; j < js[i].sprime.nPoints; j++)
       {js[i].sprime.points[j].objectives = malloc(sizeof(OBJECTIVE) * n);
        for (int k = 0; k < n; k++)
          js[i].sprime.points[j].objectives[k] = fs[0].points[j].objectives[k];
       }
     fr--;
     js[i].id = i;
     js[i].partial = inclhv(ps.points[0]);
     #if opt > 0
     qsort(js[i].sprime.points, js[i].sprime.nPoints, sizeof(POINT), greater);
     #endif
     for (int j = 0; j < js[i].sprime.nPoints; j++)
       js[i].partial -= inclhv(js[i].sprime.points[j]);
     js[i].k = js[i].sprime.nPoints - 2; 
     while (js[i].partial < 0)
       {makeDominatedBit(js[i].sprime, js[i].k);
        js[i].partial += hv(fs[0]);
        fr--;
        js[i].k--;
       }
     printf("j(%d) = %d\n", i, js[i].k);
    }
  qsort(js, ps.nPoints, sizeof(JOB), smallersofar);
  for (int i = 0; i < ps.nPoints - 1; i++) {js[i].left = i + 1; js[i].right = ps.nPoints;}
  js[ps.nPoints - 1].left = ps.nPoints; js[ps.nPoints - 1].right = ps.nPoints;
  js[ps.nPoints].partial = 99999999;
  while (js[0].k >= 0)
    {int next = js[js[0].left].partial < js[js[0].right].partial ? js[0].left : js[0].right;
     while (js[0].k >= 0 && js[0].partial <= js[next].partial)
       {makeDominatedBit(js[0].sprime, js[0].k);
        js[0].partial += hv(fs[0]);
        fr--;
        js[0].k--;
       }
     printf("j(0) = %d\n", js[0].k);
     if (js[0].partial > js[next].partial) {fixHeap(js, next); printf("loop!\n");}
    }
  return js[0].id;
}


double hv(FRONT ps)
// returns the hypervolume of ps[0 ..] 
{
  #if opt > 0
  qsort(ps.points, ps.nPoints, sizeof(POINT), greater);
  #endif

  #if opt == 2
  if (n == 2) return hv2(ps);
  #elif opt == 4
  if (n == 3) return hv3_AVL(ps);
  #else
  if (n == 3) return hv3(ps);
  #endif

  double volume = 0;

  #if opt <= 1
  for (int i = 0; i < ps.nPoints; i++) volume += exclhv(ps, i);
  #else
  n--;
  for (int i = ps.nPoints - 1; i >= 0; i--)
    // we can ditch dominated points here, 
    // but they will be ditched anyway in dominatedBit 
    volume += fabs(ps.points[i].objectives[n] - ref.objectives[n]) * exclhv(ps, i);

  n++; 
  #endif

  return volume;
}


int main(int argc, char *argv[]) 
// processes each front from the file 
{
  FILECONTENTS *f = readFile(argv[1]);

  // find the biggest fronts
  for (int i = 0; i < f->nFronts; i++)
    {if (f->fronts[i].nPoints > maxm) maxm = f->fronts[i].nPoints;
     if (f->fronts[i].n       > maxn) maxn = f->fronts[i].n;
    }

  // allocate memory
  #if opt == 0
  fs = malloc(sizeof(FRONT) * maxm);
  #else

  int optn = opt;
  // use 3 for tree 3D
  if (optn == 4)
    optn = 3;

  // slicing (opt > 1) saves a level of recursion, but in incremental mode more levels are needed
  int maxd = maxn - (optn / 2 + 1) + (optn / 3 + 1) * (1 - mode2); 
  fs = malloc(sizeof(FRONT) * maxd);

  // 3D base (opt = 3) needs space for the sentinels
  int maxp = maxm + 2 * (optn / 3);
  //int maxp = 100000;
  for (int i = 0; i < maxd; i++) 
    {fs[i].points = malloc(sizeof(POINT) * maxp); 
     for (int j = 0; j < maxp; j++) 
     {
       fs[i].points[j].tnode = malloc(sizeof(avl_node_t));
       // slicing (opt > 1) saves one extra objective at each level, but again incremental needs more 
       fs[i].points[j].objectives = malloc(sizeof(OBJECTIVE) * (maxn - (i + mode2) * (optn / 2)));
     }
    }
  #endif

    tree  = avl_alloc_tree ((avl_compare_t) compare_tree_asc,
                            (avl_freeitem_t) free);


  // initialise the reference point
  ref.objectives = malloc(sizeof(OBJECTIVE) * maxn);
  ref.tnode = malloc(sizeof(avl_node_t));
  if (argc == 2)
    {printf("No reference point provided: using the origin\n");
     for (int i = 0; i < maxn; i++) ref.objectives[i] = REFERENCE_POINT;
    }
  else if (argc - 2 != maxn)
    {printf("Your reference point should have %d values\n", maxn);
     return 0;
    }
  else 
  for (int i = 2; i < argc; i++) ref.objectives[i - 2] = atof(argv[i]);

  POINT t;
  for (int i = 0; i < f->nFronts; i++) 
    {
      
      struct timeval tv1, tv2;
      struct rusage ru_before, ru_after;
      getrusage (RUSAGE_SELF, &ru_before);

      
      n = f->fronts[i].n;
     switch (mode)
       {case 0: 
                #if opt >= 3
                if (n == 2)
                  {qsort(f->fronts[i].points, f->fronts[i].nPoints, sizeof(POINT), greater);
                   printf("hv(%d) = %1.10f\n", i+1, hv2(f->fronts[i])); 
                  }
                else
                #endif
                printf("hv(%d) = %1.10f\n",     i+1, hv(f->fronts[i])); 
                break; 
        case 1: for (int j = 0; j < f->fronts[i].nPoints; j++)
		  {t = f->fronts[i].points[0];
                   f->fronts[i].points[0] = f->fronts[i].points[j];
		   f->fronts[i].points[j] = t;
                   printf("exclhv(%d)[%d] = %1.10f\n", i+1, j+1, exclhv(f->fronts[i], 0)); 
		  }
                break; 
        case 2: printf("smallest(%d) = %d\n",   i+1, smallestpoint0(f->fronts[i])); break;
        case 3: printf("smallest(%d) = %d\n",   i+1, smallestpoint1(f->fronts[i])); break;
        case 4: printf("smallest(%d) = %d\n",   i+1, smallestpoint2(f->fronts[i])); break;
       }
          getrusage (RUSAGE_SELF, &ru_after);
          tv1 = ru_before.ru_utime;
          tv2 = ru_after.ru_utime;
      
          printf("Time: %f (s)\n", tv2.tv_sec + tv2.tv_usec * 1e-6 - tv1.tv_sec - tv1.tv_usec * 1e-6);

        //printf("Average time = %f s\n", (tv2.tv_sec + tv2.tv_usec * 1e-6 - tv1.tv_sec - tv1.tv_usec * 1e-6) / f->nFronts);

    }

  #if opt == 0
  printf("frmax = %d\n", frmax);
  #else
  printf("frmax = %d\n", maxd);
  #endif

  return 0;
}
