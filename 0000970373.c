/**
 * Elvis Perlika
 * 0000970373
 * Gruppo A
 * elvis.perlika@studio.unibo.it
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "math.h"

#define NODE_UNDEF -1
#define WET 0
#define DRY -1

typedef struct {
    int key;
    double prio;
} HeapElem;

typedef struct {
    HeapElem *heap;
    int *pos; 
    int n; 
    int size; 
} MinHeap;

/**
 * Crea un nuovo min-heap vuoto che può contenere al massimo `size` elementi.
*/
MinHeap *minheap_create(int size);

/* Svuota lo heap */
void minheap_clear(MinHeap *h);

/* Dealloca la memoria occupata dallo heap h e dal suo contenuto */
void minheap_destroy(MinHeap *h);

/* Restituisce 1 se e solo se lo heap è vuoto */
int minheap_is_empty(const MinHeap *h);

/* Restituisce 1 se e solo se lo heap è pieno */
int minheap_is_full(const MinHeap *h);

/* Ritorna il numero di elementi effettivamente presenti nello heap */
int minheap_get_n(const MinHeap *h);

/* Restituisce la chiave associata alla minima priorità; non modifica
   lo heap */
int minheap_min(const MinHeap *h);

/* Restituisce la coppia (chiave, prio) con priorità minima; non
   modifica lo heap.

   Precondizione: lo heap non deve essere vuoto. */
HeapElem minheap_min2(const MinHeap *h);

/* Inserisce una nuova chiave `key` con priorità `prio`.

   Precondizioni:
   - `key` deve essere una chiave valida;
   - `key` non deve essere già presente nello heap;
   - Lo heap non deve essere pieno. */
void minheap_insert(MinHeap *h, int key, double prio);

/* Rimuove dallo heap la coppia (chiave, prio) con priorità minima, e
   restituisce la chiave di tale coppia.

   Precondizione: lo heap non deve essere vuoto. */
int minheap_delete_min(MinHeap *h);

/* Modifica la priorità associata alla chiave `key`.

   Precondizione: la chiave `key` deve essere presente nello heap. */
void minheap_change_prio(MinHeap *h, int key, double new_prio);

/* struttura arco */
typedef struct Edge {
    int src;            /* nodo sorgente        */
    int dst;            /* nodo destinazione    */
    double weight;      /* peso dell'arco       */
    struct Edge *next;
} Edge;

typedef enum { GRAPH_UNDIRECTED, GRAPH_DIRECTED } Graph_type;

/* struttura grafo */
typedef struct {
    int n;              /* numero di nodi               */
    int m;              /* numero di archi              */
    Graph_type t;       /* tipo di grafo (orientato/non orientato) */
    Edge **edges;       /* array di liste di adiacenza  */
    int *in_deg;        /* grado entrante dei nodi      */
    int *out_deg;       /* grado uscente dei nodi       */
} Graph;

/**
 * Crea un nuovo grafo con `n` nodi e tipo `t`.
*/
Graph *graph_create(int n, Graph_type t);

/* Libera tutta la memoria occupata dal grafo e dalle liste di
   adiacenza */
void graph_destroy(Graph *g);

/* Restituisce il tipo di grafo */
Graph_type graph_type(const Graph *g);

/* Aggiunge un nuovo arco (src, dst) con peso "weight". */
void graph_add_edge(Graph *g, int src, int dst, double weight);

/* Restituisce un puntatore al primo arco della lista di adiacenza
   associata al nodo `v` (`NULL` se la lista è vuota) */
Edge *graph_adj(const Graph *g, int v);

/* Restituisce il numero di nodi del grafo */
int graph_n_nodes(const Graph *g);

/* Restituisce il numero di archi del grafo */
int graph_n_edges(const Graph *g);

Graph *graph_create( int n, Graph_type t )
{
    int i;
    Graph *g = (Graph*)malloc(sizeof(*g));
    assert(g != NULL);
    assert(n > 0);

    g->n = n;
    g->m = 0;
    g->t = t;
    g->edges = (Edge**)malloc(n * sizeof(Edge*));
    assert(g->edges != NULL);
    g->in_deg = (int*)malloc(n * sizeof(*(g->in_deg)));
    assert(g->in_deg != NULL);
    g->out_deg = (int*)malloc(n * sizeof(*(g->out_deg)));
    assert(g->out_deg != NULL);
    for (i=0; i<n; i++) {
        g->edges[i] = NULL;
        g->in_deg[i] = g->out_deg[i] = 0;
    }
    return g;
}

void graph_destroy(Graph *g)
{
    int i;

    assert(g != NULL);

    for (i=0; i<g->n; i++) {
        Edge *edge = g->edges[i];
        while (edge != NULL) {
            Edge *next = edge->next;
            free(edge);
            edge = next;
        }
        g->edges[i] = NULL; 
    }
    free(g->edges);
    free(g->in_deg);
    free(g->out_deg);
    g->n = 0;
    g->edges = NULL;
    free(g);
}

Graph_type graph_type(const Graph *g)
{
    return g->t;
}

static Edge *new_edge(int src, int dst, double weight, Edge *next)
{
    Edge *edge = (Edge*)malloc(sizeof(Edge));
    assert(edge != NULL);

    edge->src = src;
    edge->dst = dst;
    edge->weight = weight;
    edge->next = next;
    return edge;
}

static int graph_adj_insert(Graph *g, int src, int dst, double weight)
{
    g->edges[src] = new_edge(src, dst, weight, g->edges[src]);
    g->in_deg[dst]++;
    g->out_deg[src]++;
    return 0;
}

void graph_add_edge(Graph *g, int src, int dst, double weight)
{
    int status = 0;

    assert(g != NULL);

    assert((src >= 0) && (src < graph_n_nodes(g)));
    assert((dst >= 0) && (dst < graph_n_nodes(g)));

    status = graph_adj_insert(g, src, dst, weight);
    if (graph_type(g) == GRAPH_UNDIRECTED) {
        status |= graph_adj_insert(g, dst, src, weight);
    }
    if (status == 0)
        g->m++;
    else
        fprintf(stderr, "Ho ignorato l'arco duplicato (%d,%d)\n", src, dst);
}

int graph_n_nodes(const Graph *g)
{
    assert(g != NULL);

    return g->n;
}

int graph_n_edges(const Graph *g)
{
    assert(g != NULL);

    return g->m;
}

Edge *graph_adj(const Graph *g, int v)
{
    assert(g != NULL);
    assert((v >= 0) && (v < graph_n_nodes(g)));

    return g->edges[v];
}

void minheap_clear( MinHeap *h )
{
    int i;
    assert(h != NULL);
    for (i=0; i<h->size; i++) {
        h->pos[i] = -1;
    }
    h->n = 0;
}

MinHeap *minheap_create(int size)
{
    MinHeap *h = (MinHeap*)malloc(sizeof(*h));
    assert(h != NULL);
    assert(size > 0);

    h->size = size;
    h->heap = (HeapElem*)malloc(size * sizeof(*(h->heap)));
    assert(h->heap != NULL);
    h->pos = (int*)malloc(size * sizeof(*(h->pos)));
    assert(h->pos != NULL);
    minheap_clear(h);
    return h;
}

void minheap_destroy( MinHeap *h )
{
    assert(h != NULL);

    h->n = h->size = 0;
    free(h->heap);
    free(h->pos);
    free(h);
}

/* Restituisce 1 se l'indice `i` appartiene all'intervallo degli
   indici validi degli elementi validi nell'array che rappresenta lo
   heap. */
static int valid(const MinHeap *h, int i)
{
    assert(h != NULL);

    return ((i >= 0) && (i < h->n));
}

/* Scambia heap[i] con heap[j] */
static void swap(MinHeap *h, int i, int j)
{
    HeapElem tmp;

    assert(h != NULL);
    assert(valid(h, i));
    assert(valid(h, j));
    assert(h->pos[h->heap[i].key] == i);
    assert(h->pos[h->heap[j].key] == j);

    tmp = h->heap[i];
    h->heap[i] = h->heap[j];
    h->heap[j] = tmp;

    h->pos[h->heap[i].key] = i;
    h->pos[h->heap[j].key] = j;
}

/* Restituisce l'indice del padre del nodo i */
static int parent(const MinHeap *h, int i)
{
    assert(valid(h, i));

    return (i+1)/2 - 1;
}

/* Restituisce l'indice del figlio sinistro del nodo `i`. Ritorna un
   indice non valido se `i` non ha figlio sinistro. */
static int lchild(const MinHeap *h, int i)
{
    assert(valid(h, i));

    return 2*i + 1;
}

/* Restituisce l'indice del figlio destro del nodo `i`. Ritorna un
   indice non valido se `i` non ha figlio destro. */
static int rchild(const MinHeap *h, int i)
{
    assert(valid(h, i));

    return 2*i + 2;
}

/* Restituisce l'indice del figlio di `i` con priorità minima. Se `i`
   non ha figli, restituisce -1 */
static int min_child(const MinHeap *h, int i)
{
    int l, r, result = -1;

    assert(valid(h, i));

    l = lchild(h, i);
    r = rchild(h, i);
    if (valid(h, l)) {
        result = l;
        if (valid(h, r) && (h->heap[r].prio < h->heap[l].prio)) {
            result = r;
        }
    }
    return result;
}

/* Scambia l'elemento in posizione `i` con il padre fino a quando
   raggiunge la posizione corretta nello heap */
static void move_up(MinHeap *h, int i)
{
    int p;

    assert(valid(h, i));

    p = parent(h, i);
    while ( valid(h, p) && (h->heap[i].prio < h->heap[p].prio) ) {
        swap(h, i, p);
        i = p;
        p = parent(h, i);
    }
}

/* Scambia l'elemento in posizione `i` con il figlio avente priorità
   minima, fino a quando l'elemento raggiunge la posizione
   corretta. Questa funzione corrisponde a Min-Heapify() */
static void move_down(MinHeap *h, int i)
{
    int done = 0;

    assert(valid(h, i));

    do {
        const int dst = min_child(h, i);
        if (valid(h, dst) && (h->heap[dst].prio < h->heap[i].prio)) {
            swap(h, i, dst);
            i = dst;
        } else {
            done = 1;
        }
    } while (!done);
}

int minheap_is_empty(const MinHeap *h)
{
    assert(h != NULL);

    return (h->n == 0);
}

int minheap_is_full(const MinHeap *h)
{
    assert(h != NULL);

    return (h->n == h->size);
}

int minheap_get_n(const MinHeap *h)
{
    assert(h != NULL);

    return h->n;
}

int minheap_min(const MinHeap *h)
{
    assert( !minheap_is_empty(h) );

    return h->heap[0].key;
}

HeapElem minheap_min2( const MinHeap *h)
{
    assert( !minheap_is_empty(h) );

    return h->heap[0];
}

void minheap_insert(MinHeap *h, int key, double prio)
{
    int i;

    assert( !minheap_is_full(h) );
    assert((key >= 0) && (key < h->size));
    assert(h->pos[key] == -1);

    i = h->n++;
    h->pos[key] = i;
    h->heap[i].key = key;
    h->heap[i].prio = prio;
    move_up(h, i);
}

int minheap_delete_min(MinHeap *h)
{
    int result;

    assert( !minheap_is_empty(h) );

    result = minheap_min(h);
    swap(h, 0, h->n-1);
    assert( h->heap[h->n - 1].key == result );
    h->pos[result] = -1;
    h->n--;
    if (!minheap_is_empty(h)) {
        move_down(h, 0);
    }
    return result;
}

void minheap_change_prio(MinHeap *h, int key, double newprio)
{
    int j;
    double oldprio;

    assert(h != NULL);
    assert(key >= 0 && key < h->size);
    j = h->pos[key];
    assert( valid(h, j) );
    oldprio = h->heap[j].prio;
    h->heap[j].prio = newprio;
    if (newprio > oldprio) {
        move_down(h, j);
    } else {
        move_up(h, j);
    }
}

/**
 * Getter of the array index of the flattened matrix.
*/
static int getIndex(int i, int j, int nCol) {
    return i * nCol + j;
}

/* crea il metodo contrario a getIndex */
static void getIJ(int index, int nCol, int *i, int *j) {
    *i = index / nCol;
    *j = index % nCol;
}

/* funzione che legge l'array dei predecessori e ritorna il cammino dal nodo 'src' al nodo 'dst'*/
static int get_path(const int *p, int src, int dst, int *path)
{
    int i = 0;
    if (dst != NODE_UNDEF) {
        i = get_path(p, src, p[dst], path);
        path[i] = dst;
        i++;
    }
    return i;
}

/**
 * Dijkstra algorithm.
*/
void dijkstra( const Graph *g, int s, double *d, int *p)
{
    int i, u;
    Edge *adj;
    MinHeap *h = minheap_create(graph_n_nodes(g));
    for (i = 0; i < g->n; i++) {
        d[i] = HUGE_VAL;
        p[i] = NODE_UNDEF;
    }
    for (i = 0; i < g->n; i++) {
        minheap_insert(h, i, d[i]);
    }
    d[s] = 0;
    
    while (!minheap_is_empty(h))
    {
        u = minheap_delete_min(h);
        adj = graph_adj(g, u);
        while (adj != NULL)
        {
            if (d[adj->dst] > (d[u] + adj->weight)) {
                d[adj->dst] = (d[u] + adj->weight);
                p[adj->dst] = u;
                minheap_change_prio(h, adj->dst, d[adj->dst]);
            }
            adj = adj->next;
        }
    }
    minheap_destroy(h);
}

void printPath(int k, int *path, int m, int **matrix) {
    int i = 0;
    int i2 = 0; 
    int j2 = 0; 
    int w = 0;
    if(k == 1) {
        printf("-1 -1\n");
    } else {
        /* get the number of wet cells */
        w = 0;
        for (i = 0; i < k; i++) {
            getIJ(path[i], m, &i2, &j2);
            if (matrix[i2][j2] == 0) {
                w++;
            }  
        }

        /* print the path */
        printf("%d %d\n", k, w);
        for (i = 0; i < k && i + 1 < k; i++) {
            if (path[i] == path[i+1] - 1) {
                printf("E");
            } else if (path[i] == path[i+1] - m) {
                printf("S");
            } else if (path[i] == path[i+1] + 1) {
                printf("W");
            } else if (path[i] == path[i+1] + m) {
                printf("N");
            }
        }
    }
}


int **fillMatrix(FILE *filein, int n, int m) {
    int i, j, val;
    char cval;
    int **matrix;
    /* create the matrix */
    matrix = (int **)malloc(n * sizeof(int *));
    for (i = 0; i < n; i++) {
        matrix[i] = (int *)malloc(m * sizeof(int));
    }
    /* load data in the matrix */
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            fscanf(filein, "%c", &cval);
            if (cval == '\n') {
                fscanf(filein, "%c", &cval);
            }
            val = cval - '0';
            matrix[i][j] = val;
        }
    }
    return matrix;
}

void fillGraph(Graph *g, int **matrix, int n, int m) {

    int i = 0, j = 0, val = 0;

    /* set the edges with weight */
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            val = matrix[i][j];
            
            /* if the value is a sidewalk 
             * and it isn't the bot-right element of the matrix 
             * and the val's right sidewalk or val's bottom sidewalk are in the array */
            if (!(i == n - 1 && j == m - 1) && (val == WET || val == DRY)) {
                if (i + 1 < n && j == m - 1 && matrix[i+1][j] == WET) {
                    graph_add_edge(g, getIndex(i, j, m), getIndex(i+1, j, m), 1.1);
                } 
                else if (i + 1 < n && j == m - 1 && matrix[i+1][j] == DRY) {
                    graph_add_edge(g, getIndex(i, j, m), getIndex(i+1, j, m), 1);
                } 
                else if (j + 1 < m && i == n - 1 && matrix[i][j+1] == WET) {
                    graph_add_edge(g, getIndex(i, j, m), getIndex(i, j+1, m), 1.1);
                } 
                else if (j + 1 < m && i == n - 1 && matrix[i][j+1] == DRY) {
                    graph_add_edge(g, getIndex(i, j, m), getIndex(i, j+1, m), 1);
                } 
                else {
                    if (j + 1 < m && matrix[i][j+1] == WET) {
                        graph_add_edge(g, getIndex(i, j, m), getIndex(i, j+1, m), 1.1);
                    } else if (j + 1 < m && matrix[i][j+1] == DRY) {
                        graph_add_edge(g, getIndex(i, j, m), getIndex(i, j+1, m), 1);
                    }
                    if (i + 1 < n && matrix[i+1][j] == WET) {
                        graph_add_edge(g, getIndex(i, j, m), getIndex(i+1, j, m), 1.1);
                    } else if (i + 1 < n && matrix[i+1][j] == DRY) {
                        graph_add_edge(g, getIndex(i, j, m), getIndex(i+1, j, m), 1);
                    }
                }
            }
        }
    }
}

int main(int argc, char const *argv[])
{
    Graph *g = NULL;
    int n = 0, m = 0;
    int i = 0, j = 0, k = 0;
    int val = 0;
    int **matrix;
    double *d = NULL; 
    int *p = NULL;
    int *path = NULL;
    FILE *filein = stdin;
    if (argc < 2) {
        fprintf(stderr, "Usage: %s filename [src [dst]]\n", argv[0]);
        return EXIT_FAILURE;
    }

    if (strcmp(argv[1], "-") != 0) {
        filein = fopen(argv[1], "r");
        if (filein == NULL) {
            fprintf(stderr, "Can not open %s\n", argv[1]);
            return EXIT_FAILURE;
        }
    }

    fscanf(filein, "%d %d", &n, &m);
    
    /* create the graph */
    g = graph_create(n * m, GRAPH_UNDIRECTED);

    matrix = fillMatrix(filein, n, m);

    /* set the dry cells -1 */
    k = 0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            val = matrix[i][j];
            if (val > 0) {
                for (k = 1; j + k < m && k <= val; k++) {
                    if (matrix[i][j+k] == 0) {
                        matrix[i][j+k] = -1;        
                    }
                }
            } else {
                matrix[i][j] = val;
            }
        }
    }


    d = (double*)malloc((n*m) * sizeof(double));
    p = (int*)malloc((n*m) * sizeof(int));

    fillGraph(g, matrix, n, m);
    

    /* launch dijkstra */
    dijkstra(g, 0, d, p);

    /* get the path */
    path = (int*)malloc((n*m) * sizeof(int));
    /* get the size of the path */
    k = get_path(p, 0, (n*m) - 1, path);


    printPath(k, path, m, matrix);
    

    /* free memory */
    free(d);
    free(p);
    free(path);
    for (i = 0; i < n; i++) {
        free(matrix[i]);
    }
    free(matrix);
    graph_destroy(g);
    fclose(filein);

    return 0;
}


