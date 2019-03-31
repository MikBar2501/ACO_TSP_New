using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ACO_TSP_New
{
    class Parameters
    {
        public double rho = 0.98;
        public UInt32 antsCount = 10;
        public double beta = 2;
        public UInt32 candListSize = 15;
        public double pBest = 0.05;

        public double GetEvaporationRate()
        {
            return 1 - rho;
        }

        public void SetAntsCount(UInt32 ants)
        {
            antsCount = ants;
        }
    }
}
