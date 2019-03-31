using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ACO_TSP_New
{
    class Program
    {
        const UInt32 MaxCandListSize = 64;
        static void Main(string[] args)
        {
            Point[] points = TSPFileReader.ReadTspFile(@"TSP\kroA100.tsp");
            Graph graph = new Graph(points);
            Parameters parameters = new Parameters();
            parameters.SetAntsCount(graph.dimension);
            RunMMAS(graph, parameters, 1000000 / parameters.antsCount);
            Console.ReadKey();
        }

        public static Ant RunMMAS(Graph graph, Parameters parameters, UInt32 iterations)
        {
            var greedySolution = CreateSolutionNN(graph);
            var greedyCost = graph.CalculateRouteLength(greedySolution.visited);

            var initialLimits = CalcTrailLimitsMMAS(parameters, graph.dimension, greedyCost);
            var minPheromone = initialLimits.Item1;
            var maxPheromone = initialLimits.Item2;

            PheromoneMemory pheromone = new PheromoneMemory(graph.dimension, maxPheromone);
            double[] heuristic = new double[graph.dimension * graph.dimension];

            foreach(var distance in graph.edgeWeight)
           

            for(int i = 0; i < graph.dimension; i++)
            {
                for(int j = 0; j < graph.dimension; j++)
                {
                    heuristic[j + (graph.dimension * i)] = (1 / Math.Pow(graph.edgeWeight[i, j], parameters.beta));
                }
            }

            Ant[] ants = new Ant[parameters.antsCount];
            for(int i = 0; i < parameters.antsCount; i++)
            {
                ants[i] = new Ant();
            }
            Ant bestAnt = new Ant();

            for(UInt32 iteration = 0; iteration < iterations; ++iteration)
            {
                foreach(Ant ant in ants)
                {
                    ant.Initialize(graph.dimension);
                    var startNode = GetRandomUInt32(0, graph.dimension - 1);
                    ant.Visit(startNode);
                    for(UInt32 j = 1; j < graph.dimension; ++j)
                    {
                        moveAntMMAS(graph, pheromone, heuristic, ant);
                    }
                    ant.cost = graph.CalculateRouteLength(ant.visited);
                }

                var iterationBest = ants[0];
                bool newBestFound = false;

                foreach(Ant ant in ants)
                {
                    if(ant.cost < bestAnt.cost)
                    {
                        bestAnt = ant;
                        newBestFound = true;

                        Console.WriteLine("New best solution found with the cost: {0} at iteration {1}", bestAnt.cost, iteration);

                    }

                    if(ant.cost < iterationBest.cost)
                    {
                        iterationBest = ant;
                    }
                }

                if(newBestFound)
                {
                    var limits = CalcTrailLimitsMMAS(parameters, graph.dimension, bestAnt.cost);
                    minPheromone = limits.Item1;
                    maxPheromone = limits.Item2;
                }

                pheromone.EvaporateFromAll(parameters.GetEvaporationRate(), minPheromone);

                var updateAnt = iterationBest;
                double deposite = 1.0 / updateAnt.cost;
                var previousNode = updateAnt.visited[updateAnt.visited.Count - 1];
                foreach(var node in updateAnt.visited)
                {
                    pheromone.Increase(previousNode, node, deposite, maxPheromone, graph.isSymetric);
                    previousNode = node;
                }
            }
            return bestAnt;
        }

        public static double GetRandomDouble()
        {
            Random r = new Random();
            double rand = r.NextDouble();
            return rand;
        }

        public static UInt32 GetRandomUInt32(int from, UInt32 to)
        {
            Random r = new Random();
            UInt32 rand = (UInt32)r.Next(from,(int)to);
            return rand;
        }

        public static Ant CreateSolutionNN(Graph graph, UInt32 startNode = 0)
        {
            Ant ant = new Ant();
            ant.Initialize(graph.dimension);
            UInt32 currentNode = startNode;
            ant.Visit(currentNode);

            for(int i = 1; i < graph.dimension; i++)
            {
                UInt32 nextNode = NextNode(graph,(UInt32)i,ant);
                if(nextNode == currentNode)
                {
                    double minDistance = double.MaxValue;
                    for(UInt32 node = 0; node < graph.dimension; ++node)
                    {
                        if(!ant.IsVisited(node))
                        {
                            var distance = graph.edgeWeight[currentNode, node];
                            if(distance < minDistance)
                            {
                                minDistance = distance;
                                nextNode = node;
                            }
                        }
                    }
                }

                ant.Visit(nextNode);
                currentNode = nextNode;
            }
            return ant;
        }

        public static UInt32 NextNode(Graph graph, UInt32 currentNode, Ant ant)
        {

            UInt32 nextNode = currentNode;
            if(!ant.AllVisited())
            {
                double[] myNeighbours = new double[graph.dimension - 1];
                for(int i = 0; i < graph.dimension - 1;i++)
                {
                    if(i == currentNode)
                    {
                        continue;
                    }
                    myNeighbours[i] = graph.edgeWeight[(int)currentNode, i];
                }

                double min = myNeighbours.Min();

                for(int j = 0; j < graph.dimension; j++)
                {
                    if(graph.edgeWeight[currentNode,j] != min)
                    {
                        continue;
                    } else
                    {
                        nextNode = (UInt32)j;
                    }
                }
            }
            return nextNode;
        }

        public static UInt32 moveAntMMAS(Graph graph, PheromoneMemory pheromone, double [] heuristic, Ant ant)
        {
            var dimension = graph.dimension;
            var currentNode = ant.visited[ant.visited.Count - 1];
            UInt32 offset = currentNode * dimension;

            UInt32 [] candList = new UInt32 [MaxCandListSize];
            UInt32 candListSize = 0;
            for(int i = 0; i < graph.dimension; i++)
            {
                if(!ant.IsVisited((UInt32)i))
                {
                    candList[candListSize] = (UInt32)i;
                    ++candListSize;
                }
            }

            UInt32 chosenNode = currentNode;

            if(candListSize > 0)
            {
                double[] productsPrefixSum = new double[MaxCandListSize];
                double total = 0;
                for(UInt32 j = 0; j < candListSize; ++j)
                {
                    var node = candList[j];
                    var product = pheromone.Get(currentNode, node) * heuristic[offset + node];
                    total += product;
                    productsPrefixSum[j] = total;
                }

                chosenNode = candList[candListSize - 1];

                var r = GetRandomDouble() * total;
                for (UInt32 i = 0; i < candListSize; ++i)
                {
                    if(r < productsPrefixSum[i])
                    {
                        chosenNode = candList[i];
                        break;
                    }
                }
            } 
            else
            {
                double maxProduct = 0;
                for(UInt32 node = 0;node < dimension; ++node)
                {
                    if(!ant.IsVisited(node))
                    {
                        var product = pheromone.Get(currentNode, node) * heuristic[offset + node];
                        if(product > maxProduct)
                        {
                            maxProduct = product;
                            chosenNode = node;
                        }
                    }
                }
            }

            ant.Visit(chosenNode);
            return chosenNode;

        }


        public static Tuple<double,double> CalcTrailLimitsMMAS(Parameters parameters, UInt32 dimension, double solutionCost)
        {
            var tauMax = 1 / (solutionCost * (1.0 - parameters.rho));
            var avg = dimension / 2.0;
            var p = Math.Pow(parameters.pBest, 1.0 / dimension);
            var tauMin = Math.Min(tauMax, tauMax * (1 - p) / ((avg - 1) * p));
            return Tuple.Create(tauMin, tauMax);
        }


    }
}
