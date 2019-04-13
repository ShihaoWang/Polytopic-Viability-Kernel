#include "HJB_Stability.h"

// This function computes the failure metric with the value iteration method.

const double PI = 3.141592653589793238463;

// The lower and upper bound of the pendulum variables
double LLow = 0.35;             double LUpp = 1.05;
double LdotLow = -3.0;          double LdotUpp = 3.0;
double ThetaLow = -PI/2.0;      double ThetaUpp = -1.0 * ThetaLow;
double ThetadotLow = -3.0;      double ThetadotUpp = -1.0 * ThetadotLow;

const int L_Grids = 61;
const int Ldot_Grids = 61;
const int Theta_Grids = 61;
const int Thetadot_Grids = 61;

double g = 9.81;

float Ctv = 10.0;     // coefficient of theta violation
float Clv = 10.0;     // coefficient of length violation

double Clow = 0.5;
double Cupp = 2.0;

const int F_Grids = 51;       // This is the discretization of the control points within a range
double delta_t = 0.05;

// The corresponding dimension size of the StateMatrix
double L_length = LUpp - LLow;
double Ldot_length = LdotUpp - LdotLow;
double Theta_length = ThetaUpp - ThetaLow;
double Thetadot_length = ThetadotUpp - ThetadotLow;

double L_unit = L_length/(1.0* L_Grids - 1.0);
double Ldot_unit = Ldot_length/(1.0* Ldot_Grids - 1.0);
double Theta_unit = Theta_length/(1.0* Theta_Grids - 1.0);
double Thetadot_unit = Thetadot_length/(1.0* Thetadot_Grids - 1.0);

std::vector<double> L_vector(L_Grids), Ldot_vector(Ldot_Grids), Theta_vector(Theta_Grids), Thetadot_vector(Thetadot_Grids);

void StateVectorSubs(std::vector<double>& L_vector,std::vector<double>& Ldot_vector,std::vector<double>& Theta_vector, std::vector<double>& Thetadot_vector)
{
  // L_vector filling
  for (int i = 0; i < L_Grids; i++)
  {
    L_vector[i] = LLow + (1.0 * i) * L_unit;
  }
  // Ldot_vector filling
  for (int i = 0; i < Ldot_Grids; i++)
  {
    Ldot_vector[i] = LdotLow + (1.0 * i) * Ldot_unit;
  }
  // Theta_vector filling
  for (int i = 0; i < Theta_Grids; i++)
  {
    Theta_vector[i] = ThetaLow + (1.0 * i) * Theta_unit;
  }
  // Thetadot_vector filling
  for (int i = 0; i < Thetadot_Grids; i++)
  {
    Thetadot_vector[i] = ThetadotLow + (1.0 * i) * Thetadot_unit;
  }
}

SystemIndex SystemState2Index(const SystemState & x)
{
  // This function is used to transform the System State from double into the corresponding index
  // Here the x must have already been checked with the given pendulum length bound.

  double L_FloatIndex =         (x.L - LLow)/L_unit * 1.0;
  double Ldot_FloatIndex =      (x.Ldot - LdotLow)/Ldot_unit * 1.0;
  double Theta_FloatIndex =     (x.Theta - ThetaLow)/Theta_unit * 1.0;
  double Thetadot_FloatIndex =  (x.Thetadot - ThetadotLow)/Thetadot_unit * 1.0;

  int L_Index =         std::round(L_FloatIndex);
  int Ldot_Index =      std::round(Ldot_FloatIndex);
  int Theta_Index =     std::round(Theta_FloatIndex);
  int Thetadot_Index =  std::round(Thetadot_FloatIndex);

  // The L Index reshift
  if(L_Index<0)
  {
    L_Index = 0;
  }
  if(L_Index>=L_Grids)
  {
    L_Index = L_Grids-1;
  }

  // The Ldot Index reshift
  if(Ldot_Index<0)
  {
    Ldot_Index = 0;
  }
  if(Ldot_Index>=Ldot_Grids)
  {
    Ldot_Index = Ldot_Grids-1;
  }

  // The Theta Index reshift
  if(Theta_Index<0)
  {
    Theta_Index = 0;
  }
  if(Theta_Index>=Theta_Grids)
  {
    Theta_Index = Theta_Grids-1;
  }

  // The Thetadot Index reshift
  if(Thetadot_Index<0)
  {
    Thetadot_Index = 0;
  }
  if(Thetadot_Index>=Thetadot_Grids)
  {
    Thetadot_Index = Thetadot_Grids-1;
  }

  SystemIndex SystemIndex_i(L_Index, Ldot_Index, Theta_Index, Thetadot_Index);
  return SystemIndex_i;
}

SystemIndex SystemState2StateIndex(const SystemState & x)
{
  // This function is used to transform the System State from double into the corresponding index
  // Here the x must have already been checked with the given pendulum length bound.

  double L_FloatIndex =         (x.L - LLow)/L_unit * 1.0;
  double Ldot_FloatIndex =      (x.Ldot - LdotLow)/Ldot_unit * 1.0;
  double Theta_FloatIndex =     (x.Theta - ThetaLow)/Theta_unit * 1.0;
  double Thetadot_FloatIndex =  (x.Thetadot - ThetadotLow)/Thetadot_unit * 1.0;

  int L_Index =         std::round(L_FloatIndex);
  int Ldot_Index =      std::round(Ldot_FloatIndex);
  int Theta_Index =     std::round(Theta_FloatIndex);
  int Thetadot_Index =  std::round(Thetadot_FloatIndex);

  // The L Index reshift
  if(L_Index<0)
  {
    L_Index = 0;
  }
  if(L_Index>=L_Grids)
  {
    L_Index = L_Grids-1;
  }

  // The Ldot Index reshift
  if(Ldot_Index<0)
  {
    Ldot_Index = 0;
  }
  if(Ldot_Index>=Ldot_Grids)
  {
    Ldot_Index = Ldot_Grids-1;
  }

  // The Theta Index reshift
  if(Theta_Index<0)
  {
    Theta_Index = 0;
  }
  if(Theta_Index>=Theta_Grids)
  {
    Theta_Index = Theta_Grids-1;
  }

  // The Thetadot Index reshift
  if(Thetadot_Index<0)
  {
    Thetadot_Index = 0;
  }
  if(Thetadot_Index>=Thetadot_Grids)
  {
    Thetadot_Index = Thetadot_Grids-1;
  }

  SystemIndex SystemIndex_i(L_Index, Ldot_Index, Theta_Index, Thetadot_Index);
  return SystemIndex_i;
}

void FailureMetricUpdate(DBNode& Node_kp1)
{
  // KE obj
  double L_kp1, Ldot_kp1, Theta_kp1, Thetadot_kp1;
  L_kp1 = Node_kp1.NodeState.L;
  Ldot_kp1 = Node_kp1.NodeState.Ldot;
  Theta_kp1 = Node_kp1.NodeState.Theta;
  Thetadot_kp1 = Node_kp1.NodeState.Thetadot;

  double Node_k_L = Node_kp1.NodeState.L - Node_kp1.NodeState.Ldot * delta_t;


  if((Theta_kp1>=0)&&(Thetadot_kp1>0))
  {
    // On the positive side and velocity is positive
    Node_kp1.FailureMetric = 0.0;
    return;
  }

  float KE_ijkl, SC_ijkl, IC_ijkl;

  /*
    Term1: Kinetic Energy Cost
  */
  KE_ijkl = Ldot_kp1 * Ldot_kp1  + L_kp1 * L_kp1 * Thetadot_kp1 * Thetadot_kp1;

  /*
    Term2: Side Cost
  */
  switch (Node_kp1.ThetaViableFlag)
  {
    case 0:
    SC_ijkl = -1.0 * Ctv * Theta_kp1;
    break;
    case 1:
    SC_ijkl = 0.0;
    break;
    default:
    break;
  }

  /*
    Term3: Infeasibility Cost
  */
  switch (Node_kp1.LengthFeasibleFlag)
  {
    case 0:
    // Length Infeasible FailureMetric
    if(Node_k_L<LLow)
    {
      IC_ijkl = Clv * (LLow - Node_k_L);
    }
    else
    {
      IC_ijkl = Clv * (Node_k_L - LUpp);
    }
    break;
    case 1:
    IC_ijkl = 0.0;
    break;
    default:
    break;
  }

  // Total objective function
  Node_kp1.FailureMetric = KE_ijkl + SC_ijkl + IC_ijkl;
}

void DBNodeTLcUpdate(DBNode & Node_kp1)
{
  // This function can only be called after Node_kp1 has been initialized
  // The following properties for Node_kp1 will be updated.
  //  1. ThetaViableFlag
  //  2. LengthFeasibleFlag
  //  3. FailureMetric

  // 1. ThetaViableFlag
  if(Node_kp1.NodeState.Theta<0){  Node_kp1.ThetaViableFlag = 0;  }
  // 2. LengthFeasibleFlag
  double Node_k_L = Node_kp1.NodeState.L - Node_kp1.NodeState.Ldot * delta_t;
  if((Node_k_L<LLow)||(Node_k_L>LUpp))
  {
    Node_kp1.LengthFeasibleFlag = 0;
  }
  // 3. FailureMetric
  FailureMetricUpdate(Node_kp1);
}

double Backward_Evaluation(DBNode& Node_kp1, Eigen::Tensor<DBNode,4>& StateNodeMatrix)
{
  // This function conducts the first-order Euler integration from the k+1th-layer node to the kth layer node.
  /*
      L_k = L_(k+1) - Ldot_(k+1) * delta_t;
      Ldot_k = Ldot_(k+1) - Lddot_(k+1) * delta_t;
      Theta_k = Theta_(k+1) - Thetadot_(k+1) * delta_t;
      Thetadot_k = Thetadot_(k+1) - Thetaddot_(k+1) * delta_t;
  */

  std::set< std::tuple<int, int, int, int>> Reachable_Set;

  // kth layer
  double L_kp1=           Node_kp1.NodeState.L;
  double Ldot_kp1 =       Node_kp1.NodeState.Ldot;
  double Theta_kp1 =      Node_kp1.NodeState.Theta;
  double Thetadot_kp1 =   Node_kp1.NodeState.Thetadot;
  double Thetaddot_kp1 =  g/L_kp1 * sin(Theta_kp1) - 2.0 * Thetadot_kp1 * Ldot_kp1/L_kp1;

  // kth layer
  double L_k, Ldot_k, Theta_k, Thetadot_k;

  // Backward Integration
  switch (Node_kp1.LengthFeasibleFlag)
  {
    case 0:
    {
      // This means that there is no need to integrate this node
      return 0.0;
    }
    break;
    default:
    {
      L_k = L_kp1 - Ldot_kp1 * delta_t;
    }
    break;
  }

  Theta_k = Theta_kp1 - Thetadot_kp1 * delta_t;
  Thetadot_k = Thetadot_kp1 - Thetaddot_kp1 * delta_t;

  // The bounds on the acceleration of Pendulum Length
  double Lddot_Low = L_kp1 * Thetadot_kp1 * Thetadot_kp1 - g * cos(Theta_kp1);
  Lddot_Low = Clow * Lddot_Low;
  double Lddot_Upp = Cupp * g * cos(Theta_kp1)/(LLow - LUpp) * (L_kp1 - LLow) + Cupp * g * cos(Theta_kp1);

  double Lddot_Min = Lddot_Low;
  double Lddot_Max = max(Lddot_Low, Lddot_Upp);
  std::vector<double> Lddot_vec = linspace(Lddot_Min, Lddot_Max, F_Grids);

  // Ldot integration
  for (int i = 0; i < F_Grids; i++)
  {
    Ldot_k = Ldot_kp1 - Lddot_vec[i] * delta_t;     // Now Lddot_vec saves the acceleration of pendulum length at k+1th node
    SystemState x_k(L_k, Ldot_k, Theta_k, Thetadot_k);
    SystemIndex x_k_index = SystemState2StateIndex(x_k);
    Reachable_Set.insert(make_tuple(x_k_index.L_index, x_k_index.Ldot_index, x_k_index.Theta_index, x_k_index.Thetadot_index));
  }

  double FailureMetricVia = 0.0;
  switch (Reachable_Set.size())
  {
    case 0:
    {
      // This means that the failure metric for this node has not been changed due to the infeasibility.
      return FailureMetricVia;
    }
    default:
    {
      const int ReachableNodeNumber = Reachable_Set.size();
      std::vector<double> Reachable_FailureMetric(ReachableNodeNumber);
      std::set<std::tuple<int, int, int, int>>::iterator Reachable_Set_Itr = Reachable_Set.begin();
      for (int i = 0; i < ReachableNodeNumber; i++)
      {
        switch (i)
        {
          case 0:
          break;
          default:
          std::advance(Reachable_Set_Itr, 1);
          break;
        }
        // Here DBNodePtr points to the node at kth layer
        DBNode* DBNodePtr =&StateNodeMatrix(std::get<0>(*Reachable_Set_Itr),std::get<1>(*Reachable_Set_Itr),std::get<2>(*Reachable_Set_Itr),std::get<3>(*Reachable_Set_Itr));
        if(DBNodePtr->FailureMetric > Node_kp1.FailureMetric)
        {
          // Then this means that there is a better solution for the node at kth layer, so it should point to the current ndoe at k+1th layer.
          Reachable_FailureMetric[i] = DBNodePtr->FailureMetric - Node_kp1.FailureMetric;
          DBNodePtr->FailureMetric = Node_kp1.FailureMetric;
          DBNodePtr->NextNodeIndex = Node_kp1.SelfIndex;
        }
      }
      FailureMetricVia = *std::max_element(Reachable_FailureMetric.begin(), Reachable_FailureMetric.end());
    }
    break;
  }
  return FailureMetricVia;
}

void HJBDataWriter(const std::vector<int>& NodeParentIndVector, const std::vector<float>& NodeFailureMetricVector, const std::vector<double> & StateVectorSpecs)
{
  // *. Node Transition
  FILE * StateTransFile = NULL;
  StateTransFile = fopen("PVKNextInd.bin", "wb");
  fwrite(&NodeParentIndVector[0], sizeof(int), NodeParentIndVector.size(), StateTransFile);
  fclose(StateTransFile);

  // *. FailureMetric
  FILE * StateFailureMetricFile = NULL;
  StateFailureMetricFile = fopen("PVKFailureMetric.bin", "wb");
  fwrite(&NodeFailureMetricVector[0], sizeof(float), NodeFailureMetricVector.size(), StateFailureMetricFile);
  fclose(StateFailureMetricFile);


  FILE * StateVectorSpecsFile = NULL;
  StateVectorSpecsFile = fopen("PVkDataSpecs.bin", "wb");
  fwrite(&StateVectorSpecs[0], sizeof(double), StateVectorSpecs.size(), StateVectorSpecsFile);
  fclose(StateVectorSpecsFile);
}

void StateNodeMatrixInit(Eigen::Tensor<DBNode,4> &StateNodeMatrix, std::vector<float> &NodeFailureMetricVector)
{
  // This function is used to initialize the State Node Matrix
  DBNode * Node_ptr;
  for (int i = 0; i < L_Grids; i++)
  {
    for (int j = 0; j < Ldot_Grids; j++)
    {
      for (int k = 0; k < Theta_Grids; k++)
      {
        for (int l = 0; l < Thetadot_Grids; l++)
        {
          Node_ptr = &StateNodeMatrix(i,j,k,l);
          Node_ptr->InitUpdate(i, j, k, l, L_Grids, Ldot_Grids, Theta_Grids, Thetadot_Grids, L_vector, Ldot_vector, Theta_vector, Thetadot_vector);
          DBNodeTLcUpdate(*Node_ptr);
          NodeFailureMetricVector[Node_ptr->SelfIndex] = Node_ptr->FailureMetric;
        }
      }
    }
  }
}

int main()
{
  std::vector<double> StateVectorSpecs;

  StateVectorSpecs.push_back(LLow);
  StateVectorSpecs.push_back(LUpp);
  StateVectorSpecs.push_back(LdotLow);
  StateVectorSpecs.push_back(LdotUpp);
  StateVectorSpecs.push_back(ThetaLow);
  StateVectorSpecs.push_back(ThetaUpp);
  StateVectorSpecs.push_back(ThetadotLow);
  StateVectorSpecs.push_back(ThetadotUpp);

  StateVectorSpecs.push_back(L_Grids * 1.0);
  StateVectorSpecs.push_back(Ldot_Grids * 1.0);
  StateVectorSpecs.push_back(Theta_Grids * 1.0);
  StateVectorSpecs.push_back(Thetadot_Grids * 1.0);

  // LLow and LUpp are appended again for historical reasons.
  StateVectorSpecs.push_back(LLow);
  StateVectorSpecs.push_back(LUpp);
  StateVectorSpecs.push_back(g);

  /*
    Main computation objects initialization
  */
  StateVectorSubs(L_vector, Ldot_vector, Theta_vector, Thetadot_vector);
  Eigen::Tensor<DBNode,4> StateNodeMatrix(L_Grids, Ldot_Grids, Theta_Grids, Thetadot_Grids);

  std::vector<int> NodeTransIndVector(L_Grids * Ldot_Grids * Theta_Grids * Thetadot_Grids);
  std::vector<float> NodeFailureMetricVector(L_Grids * Ldot_Grids * Theta_Grids * Thetadot_Grids);
  StateNodeMatrixInit(StateNodeMatrix, NodeFailureMetricVector);

  // #TODO g: projection has to be discretized later

  std::clock_t start; double duration; start = std::clock();

  double FailureMetricVia = 0.1;              // Here this value is the evaluation of the failure metric change to determine the convergence.
  double FailureMetricTol = 1e-10;            // Tolerence for convergence completeness
  double FailureMetricVia_i;

  while(FailureMetricVia>FailureMetricTol)
  {
    // Now this inner for loop turns out to be very easy to be solved with value iteration.
    double FailureMetricVia_ref = 0.0;
    for (int i = 0; i < L_Grids; i++)
    {
      for (int j = 0; j < Ldot_Grids; j++)
      {
        for (int k = 0; k < Theta_Grids; k++)
        {
          for (int l = 0; l < Thetadot_Grids; l++)
          {
            double FailureMetricVia_i = Backward_Evaluation(StateNodeMatrix(i,j,k,l), StateNodeMatrix);
            // double FailureMetricVia_i = Backward_Evaluation(StateNodeMatrix(8,30,31,30), StateNodeMatrix);
            if(FailureMetricVia_i>FailureMetricVia_ref)
            {
              FailureMetricVia_ref = FailureMetricVia_i;
            }
          }
        }
      }
    }
    FailureMetricVia = FailureMetricVia_ref;
    printf ("FailureMetricVia value: %f \n", FailureMetricVia);
  }

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  printf ("Total running time: %f s\n", duration);

  int IterIndex = 0;
  for (int i = 0; i < L_Grids; i++)
  {
    for (int j = 0; j < Ldot_Grids; j++)
    {
      for (int k = 0; k < Theta_Grids; k++)
      {
        for (int l = 0; l < Thetadot_Grids; l++)
        {
          NodeTransIndVector[IterIndex] = StateNodeMatrix(i,j,k,l).NextNodeIndex;
          NodeFailureMetricVector[IterIndex] = StateNodeMatrix(i,j,k,l).FailureMetric;
          IterIndex = IterIndex + 1;
        }
      }
    }
  }

  HJBDataWriter(NodeTransIndVector, NodeFailureMetricVector, StateVectorSpecs);
  return 0;
}
