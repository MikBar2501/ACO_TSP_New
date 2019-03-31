using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ACO_TSP_New
{
    class PheromoneMemory
    {
        UInt32 dimension;
        double[] pheromoneValues;
        double minPheromoneValue;

        public PheromoneMemory(UInt32 dim, double minPherValues = 0)
        {
            minPheromoneValue = minPherValues;
            dimension = dim;
            pheromoneValues = new double[dimension * dimension];
            for (int i = 0; i < pheromoneValues.Length; i++)
            {
                pheromoneValues[i] = minPherValues;
            }
        }

        public double Get(UInt32 from, UInt32 to)
        {
            return pheromoneValues[from * dimension + to];
        }

        public void EvaporateFromAll(double evaporationRate, double minPheromoneValue)
        {
            for (int i = 0; i < pheromoneValues.Length; i++)
            {
                pheromoneValues[i] = Math.Max(minPheromoneValue, pheromoneValues[i] * (1 - evaporationRate));
            }
        }

        public void Increase(UInt32 from, UInt32 to, double deposite, double maxPheromoneValue, bool isSymmetric)
        {
            var pher = pheromoneValues[from * dimension + to];
            pher = Math.Min(maxPheromoneValue, pher + deposite);

            if (isSymmetric)
            {
                pheromoneValues[to * dimension + from] = pher;
            }
        }

    }
}
