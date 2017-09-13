#ifndef PYNE_DQKIQSJ4SNG7VAB5LX36BLIYMA
#define PYNE_DQKIQSJ4SNG7VAB5LX36BLIYMA

#include <map>
#include <vector>

namespace pyne {
namespace transmuters {


/// Basic CRAM solver that takes a (flat) A matrix and an inital nuclide
/// atom fraction composition map and returns the value.
/// \param A The transmutation matrix [unitless]
/// \param n0 The initial compositions [atom fraction]
/// \param order The order of approximation, default 14.
/// \return n1 The result of the transmutation [atom fraction]
std::map<int, double> cram(std::vector<double>& A,
                           const std::map<int, double>& n0,
                           const int order=14);

} // namespace transmuters
} // namespace pyne
#endif // PYNE_DQKIQSJ4SNG7VAB5LX36BLIYMA