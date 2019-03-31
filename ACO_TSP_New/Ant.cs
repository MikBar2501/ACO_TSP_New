using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ACO_TSP_New
{
    class Ant
    {
        public List<UInt32> visited;
        public bool[] isVisited;
        public double cost = double.MaxValue;

        public void Initialize(UInt32 dimension)
        {
            visited = new List<UInt32>();
            isVisited = new bool[dimension];
            for (int i = 0; i < isVisited.Length; i++)
            {
                isVisited[i] = false;
            }
        }

        public void Visit(UInt32 node)
        {
            visited.Add(node);
            isVisited[(int)node] = true;
        }

        public bool IsVisited(UInt32 node)
        {
            return isVisited[(int)node];
        }

        public bool AllVisited()
        {
            foreach (bool node in isVisited)
            {
                if (node == false)
                {
                    return false;
                }
            }
            return true;
        }
    }
}
