using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;

namespace ACO_TSP_New
{
    public struct Point
    {
        public Vector2 point;
        public int id;

        public Point(int id, Vector2 point)
        {
            this.id = id;
            this.point = point;
        }
    }
    class Graph
    {
        public UInt32 dimension;
        public bool isSymetric;
        public double[,] edgeWeight;
        public Point[] points;

        public Graph(Point[] points, bool isSymetric = true)
        {
            this.points = points;
            this.isSymetric = isSymetric;
            this.dimension = (UInt32)points.Length;
            edgeWeight = SetWeight(points);

        }

        public double[,] SetWeight(Point[] points)
        {
            double[,] weightArray = new double[points.Length, points.Length];
            for (int i = 0; i < points.Length; i++)
            {
                for (int j = 0; j < points.Length; j++)
                {
                    if (i == j)
                    {
                        weightArray[i, j] = 0;
                    }
                    else
                    {
                        weightArray[i, j] = Vector2.Distance(points[i].point, points[j].point);
                    }
                }
            }
            return weightArray;
        }

        /*public UInt32 FindNearestNeighbour(int index)
        {
            double min;
            UInt32 nearest;
            if(index == 0)
            {
                min = edgeWeight[index, 1];
                nearest = 1;
            } else
            {
                min = edgeWeight[index, 0];
                nearest = 0;
            }

            for(int i = 0; i< dimension; i++)
            {
                if(i == index)
                {
                    continue;
                }
                if(edgeWeight[index,i] < min)
                {
                    min = edgeWeight[index, i];
                    nearest = (UInt32)i;
                }
            }
            return nearest;
            
        }*/

        public double GetDistances(UInt32 from, UInt32 to)
        {
            return edgeWeight[from,to];
        }

        public double CalculateRouteLength(List<UInt32> route)
        {
            double distance = 0;
            if(!(route == null))
            {
                var previousNode = route[route.Count - 1];
                foreach(var node in route)
                {
                    distance += GetDistances(previousNode, node);
                    previousNode = node;
                }
            }
            return distance;
        }
    }
}
