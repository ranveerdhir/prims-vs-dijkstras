/* PrimVsKruskal.java
   
   
   The file includes the "import edu.princeton.cs.algs4.*;" so that you can use
   any of the code in the algs4.jar file. You should be able to compile your program
   with the command
   
	javac -cp .;algs4.jar PrimVsKruskal.java
	
   To conveniently test the algorithm with a large input, create a text file
   containing a test graphs (in the format described below) and run
   the program with
   
	java -cp .;algs4.jar PrimVsKruskal file.txt
	
   where file.txt is replaced by the name of the text file. Note: different operating systems have different commands.
   
   
   The input consists of a graph (as an adjacency matrix) in the following format:
   
    <number of vertices>
	<adjacency matrix row 1>
	...
	<adjacency matrix row n>
	
   Entry G[i][j] >= 0.0 of the adjacency matrix gives the weight (as type double) of the edge from 
   vertex i to vertex j (if G[i][j] is 0.0, then the edge does not exist).
   Note that since the graph is undirected, it is assumed that G[i][j]
   is always equal to G[j][i].
*/

 import edu.princeton.cs.algs4.*;

import java.util.Arrays;
import java.util.Scanner;
 import java.io.File;

//
public class PrimVsKruskal{

	/* PrimVsKruskal(G)
		Given an adjacency matrix for connected graph G, with no self-loops or parallel edges,
		determines if the minimum spanning tree of G found by Prim's algorithm is equal to 
		the minimum spanning tree of G found by Kruskal's algorithm.
		
		If G[i][j] == 0.0, there is no edge between vertex i and vertex j
		If G[i][j] > 0.0, there is an edge between vertices i and j, and the
		value of G[i][j] gives the weight of the edge.
		No entries of G will be negative.
	*/
	 
	static boolean PrimVsKruskal(double[][] G){
		int n = G.length;
		
		

		/* Build the MST by Prim's and the MST by Kruskal's */
	

         class PrimMST {
            private static final double FLOATING_POINT_EPSILON = 1.0E-12;

    		private int V;              // Number of vertices
    		private double[] distTo;    // distTo[v] = weight of shortest edge from tree to non-tree vertex
        	private boolean[] marked;   // marked[v] = true if v on tree, false otherwise
    		private IndexMinPQ<Double> pq;
			private int[] edgeTo;
			private UF uf;

    /**
     * Compute a minimum spanning tree (or forest) of an edge-weighted graph represented as an adjacency matrix.
     * @param G the weighted adjacency matrix
     */
   			 public PrimMST(double[][] G) {
		
        		V = G.length;
        		distTo = new double[V];
        		marked = new boolean[V];
        		pq = new IndexMinPQ<>(V);
				edgeTo = new int[V]; 
				uf = new UF(V);

        		for (int v = 0; v < V; v++)
            		distTo[v] = Double.POSITIVE_INFINITY;

        		for (int v = 0; v < V; v++)
            		if (!marked[v])
                    prim(G, v); // minimum spanning forest


    
    }

    // Run Prim's algorithm in graph G, starting from vertex s
    public void prim(double[][] G, int s) {
        distTo[s] = 0.0;
        pq.insert(s, distTo[s]);

        while (!pq.isEmpty()) {
            int v = pq.delMin();
			scan(G,v);	

        }
    }
	public void  scan(double[][] G, int v) {
        marked[v] = true;
	
		for (int w = 0; w < V; w++) {
			if (!marked[w] && G[v][w] > 0.0) {
				if (G[v][w] < distTo[w] && uf.find(v) != uf.find(w)) {
					distTo[w] = G[v][w];
					edgeTo[w] = v;
					if (pq.contains(w)) {
						pq.decreaseKey(w, distTo[w]);
					} else {
						pq.insert(w, distTo[w]);
					}
				}
			}
		}
	}
    // Helper method to update the distTo[] array
    
    public Queue<Edge> getMSTEdges() {
		Queue<Edge> mstEdges = new Queue<>();
		for (int v = 1; v < V; v++) {
			int w = edgeTo[v];
			mstEdges.enqueue(new Edge(w, v, G[w][v]));
		}
		return mstEdges;
	}
   

    // Print the Minimum Spanning Tree
    private void printMST(double[][] G) {
       double f = 0;
        for (int v = 1; v < V; v++) {
			int w = edgeTo[v];
            System.out.println(w + " - " + v+" "+ G[w][v]);
			f = f + G[w][v];
        }
		System.out.println(f);
    }

}
	 
    // Other methods remain the same
    // ...		
			// Instantiate the PrimMST class
			
		
			// Call the printMST method from PrimMST class
			

			
		
			// Rest of the code...
			 class KruskalMST {
				public double weight;  // Weight of MST
				private Queue<Edge> mst = new Queue<>(); // Edges in MST
			
				public KruskalMST(double[][] G) {
					int V = G.length;
			
					// Create edges from the adjacency matrix
					Edge[] edges = new Edge[V * (V - 1) / 2]; // Maximum number of edges in an undirected graph
					int edgeCount = 0;
					for (int i = 0; i < V; i++) {
						for (int j = i + 1; j < V; j++) {
							if (G[i][j] > 0.0) {
								edges[edgeCount++] = new Edge(i, j, G[i][j]);
							}
						}
					}
					Arrays.sort(edges, 0, edgeCount);
			
					// Kruskal's algorithm
					UF uf = new UF(V);
					for (int i = 0; i < edgeCount && mst.size() < V - 1; i++) {
						Edge e = edges[i];
						int v = e.either();
						int w = e.other(v);
			
						// v-w does not create a cycle
						if (uf.find(v)!=uf.find(w)) {
							uf.union(v, w); // Merge v and w components
							mst.enqueue(e); // Add edge e to mst
							weight += e.weight();
						}
					}
			
				
			
				}
				public Queue<Edge> getMSTEdges() {
					return mst;
				}


				public void printMST() {
					double k =0;
					for (Edge e : mst) {
						int v = e.either();
						int w = e.other(v);
						double weight = e.weight();
						System.out.println(v + " - " + w + " " + weight);
						k = k+ weight;
					}
					System.out.println(k);
				}
				
				
			
			
			}
			


			
			
			PrimMST primsMST = new PrimMST(G);
			KruskalMST kruskalsMST = new KruskalMST(G);
		
			// Retrieve MST edges from both algorithms
			Queue<Edge> mst1 = primsMST.getMSTEdges();
			Queue<Edge> mst2 = kruskalsMST.getMSTEdges();

			if (mst1.size() != mst2.size()) {
				return false;
			}
			Queue<Edge> unvisitedEdges = new Queue<Edge>();
			for (Edge edge : mst2) {
				unvisitedEdges.enqueue(edge);
			}
			for (Edge edge1 : mst1) {
				boolean found = false;
				Queue<Edge> tempQueue = new Queue<>();
	
				for (Edge edge2 : unvisitedEdges) {
					if ((edge1.either() == edge2.either() && edge1.other(edge1.either()) == edge2.other(edge2.either()))
							|| (edge1.either() == edge2.other(edge2.either()) && edge1.other(edge1.either()) == edge2.either())) {
						found = true;
					} else {
						tempQueue.enqueue(edge2);
					}
				}
	
				if (!found) {
					return false;
				}
	
				unvisitedEdges = tempQueue;
			}
			return true;
			
	
		
			// Compare edges and their weights
			
	}
		
	/* main()
	   Contains code to test the PrimVsKruskal function. You may modify the
	   testing code if needed. 
	*/
   public static void main(String[] args) {
		Scanner s;
		if (args.length > 0){
			try{
				s = new Scanner(new File(args[0]));
			} catch(java.io.FileNotFoundException e){
				System.out.printf("Unable to open %s\n",args[0]);
				return;
			}
			System.out.printf("Reading input values from %s.\n",args[0]);
		}else{
			s = new Scanner(System.in);
			System.out.printf("Reading input values from stdin.\n");
		}
		
		int n = s.nextInt();
		double[][] G = new double[n][n];
		int valuesRead = 0;
		for (int i = 0; i < n && s.hasNextDouble(); i++){
			for (int j = 0; j < n && s.hasNextDouble(); j++){
				G[i][j] = s.nextDouble();
				if (i == j && G[i][j] != 0.0) {
					System.out.printf("Adjacency matrix contains self-loops.\n");
					return;
				}
				if (G[i][j] < 0.0) {
					System.out.printf("Adjacency matrix contains negative values.\n");
					return;
				}
				if (j < i && G[i][j] != G[j][i]) {
					System.out.printf("Adjacency matrix is not symmetric.\n");
					return;
				}
				valuesRead++;
			}
		}
		
		if (valuesRead < n*n){
			System.out.printf("Adjacency matrix for the graph contains too few values.\n");
			return;
		}	
		
        boolean pvk = PrimVsKruskal(G);
        System.out.printf("Does Prim MST = Kruskal MST? %b\n", pvk);
    }
}
