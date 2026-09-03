#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>

#include "Radial_Binning.hpp"

namespace
{

constexpr double R_SUN_KM = 6.957e5;
constexpr double R_EARTH_KM = 6371.0;

TEST(BincountRadialGrid, PreservesLegacySolarDefaultGeometry)
{
    DaMaSCUS_SUN::Bincount_Radial_Grid grid(R_SUN_KM);

    EXPECT_EQ(1100u, grid.Inner_Bin_Count());
    EXPECT_EQ(423u, grid.Exterior_Bin_Count());
    EXPECT_EQ(1523u, grid.Bin_Count());

    EXPECT_NEAR(695.7, grid.Inner_Bin_Width_Km(), 1.0e-10);
    EXPECT_NEAR(1.1 * R_SUN_KM, grid.Inner_Extent_Km(), 1.0e-8);

    EXPECT_NEAR(
        DaMaSCUS_SUN::BINCOUNT_RADIAL_GRID_AU_KM,
        grid.Domain_Max_Km(),
        1.0e-6);

    EXPECT_EQ(0, grid.Bin_Index_Km(0.0));
    EXPECT_EQ(1099, grid.Bin_Index_Km(grid.Inner_Extent_Km() - 1.0e-6));
    EXPECT_EQ(1100, grid.Bin_Index_Km(grid.Inner_Extent_Km()));
    EXPECT_EQ(-1, grid.Bin_Index_Km(grid.Domain_Max_Km()));
}

TEST(BincountRadialGrid, EarthUsesEarthRadiusRatherThanSolarRadius)
{
    DaMaSCUS_SUN::Bincount_Radial_Grid grid(R_EARTH_KM);

    EXPECT_EQ(1100u, grid.Inner_Bin_Count());

    EXPECT_NEAR(6.371, grid.Inner_Bin_Width_Km(), 1.0e-12);
    EXPECT_NEAR(7008.1, grid.Inner_Extent_Km(), 1.0e-10);

    // At r = 0.5 R_earth, the 0.001 R_earth grid must select bin 500.
    EXPECT_EQ(500, grid.Bin_Index_Km(0.5 * R_EARTH_KM));

    // Earth has much smaller exterior shells than the Sun, hence it needs
    // more shells before reaching the same 1 AU outer boundary.
    EXPECT_GT(grid.Exterior_Bin_Count(), 423u);

    const std::vector<double>& edges = grid.Bin_Edges_Km();

    ASSERT_EQ(grid.Bin_Count() + 1, edges.size());
    EXPECT_DOUBLE_EQ(0.0, edges.front());
    EXPECT_DOUBLE_EQ(grid.Domain_Max_Km(), edges.back());

    for(std::size_t i = 1; i < edges.size(); i++)
        EXPECT_GT(edges[i], edges[i - 1]);
}

} // namespace
