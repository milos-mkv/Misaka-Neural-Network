#include "MisakaNeuralNetwork.hpp"

int main(const int argc, const char** argv)
{
    Misaka::NeuralNetwork neuralNet(2, 2, 1);

    // Should train multiple times in random order.
    neuralNet.Train({ 0, 0 }, { 0 });
    neuralNet.Train({ 0, 1 }, { 1 });
    neuralNet.Train({ 1, 0 }, { 1 });
    neuralNet.Train({ 1, 1 }, { 0 });
    /// ...
    /// 
    neuralNet.Query({ 0, 0 }).Print();  // 0.066057
    neuralNet.Query({ 0, 1 }).Print();  // 0.923139
    neuralNet.Query({ 1, 0 }).Print();  // 0.937463
    neuralNet.Query({ 1, 1 }).Print();  // 0.062518

    return 0;
}