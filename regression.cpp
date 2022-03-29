#include "regression.hpp"

int Regression::size() const { return points_.size();}
void Regression::add(double x, double y) {
  points_.push_back ({x,y});
  } 
Result Regression::fit() const {
  double sum_x_{};
  double sum_y_{};
  double sum_xy_{};
  double sum_x2_{};
  for (auto const& p : points_) {
    sum_x_ += p.x;
    sum_y_ += p.y;
    sum_xy_ += p.x * p.y;
    sum_x2_ += p.x * p.x;
  }
  int N_ = points_.size();

  if (N_ < 2) {
      throw std::runtime_error{"Not enough points to run a fit"};
    }
    double const d = N_ * sum_x2_ - sum_x_ * sum_x_;
    if (d == 0.) {
      throw std::runtime_error{"Trying to fit a vertical line"};
    }
    double const a = (sum_y_ * sum_x2_ - sum_x_ * sum_xy_) / d;
    double const b = (N_ * sum_xy_ - sum_x_ * sum_y_) / d;
    return {a, b};
    }
Result fit(Regression const &reg) {
  // the following call would fail compilation if the fit method weren't const
  return reg.fit(); }
