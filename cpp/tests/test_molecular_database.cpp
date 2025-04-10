#include <gtest/gtest.h>
#include "MolecularDatabase.h"
#include "ForceField.h"
#include "AtomicConfiguration.h"

class MolecularDatabaseTest : public ::testing::Test {
protected:
    MolecularDatabase db;
    ForceField mol1, mol2;
    Mat3d lat_vec;

    void SetUp() override {
        // Set up mol1 - square molecule
        mol1.natoms = 4;
        mol1.apos = new Vec3d[4];
        mol1.atypes = new int[4];
        mol1.apos[0].set(0.0, 0.0, 0.0);    // First atom at origin
        mol1.apos[1].set(1.0, 0.0, 0.0);    // Second atom 1Å away in x direction
        mol1.apos[2].set(1.0, 1.0, 0.0);    // Third atom forms corner of square
        mol1.apos[3].set(0.0, 1.0, 0.0);    // Fourth atom completes the square
        for(int i = 0; i < 4; i++) mol1.atypes[i] = 1;

        // Set up mol2 - initially same as mol1 but shifted
        mol2.natoms = 4;
        mol2.apos = new Vec3d[4];
        mol2.atypes = new int[4];
        mol2.apos[0].set(5.0, 0.0, 0.0);
        mol2.apos[1].set(6.0, 0.0, 0.0);
        mol2.apos[2].set(6.0, 1.0, 0.0);
        mol2.apos[3].set(5.0, 1.0, 0.0);
        for(int i = 0; i < 4; i++) mol2.atypes[i] = 1;

        // Set up periodic boundary conditions
        lat_vec.vecs[0].set(1.0, 0.0, 0.0);  // periodic in x with period 1
        lat_vec.vecs[1].set(0.0, 1.0, 0.0);  // periodic in y with period 1
        lat_vec.vecs[2].set(0.0, 0.0, 1.0);  // periodic in z with period 1
    }

    // void TearDown() override {
    //     delete[] mol1.apos;
    //     delete[] mol1.atypes;
    //     delete[] mol2.apos;
    //     delete[] mol2.atypes;
    // }
};

// Test 1: Identical molecules with periodic boundary conditions
TEST_F(MolecularDatabaseTest, IdenticalMolecules) {
    double distance = db.computeDistanceOnSurf(&mol1, &mol2, &lat_vec);
    EXPECT_LT(distance, 0.1) << "Distance between identical molecules should be near zero with periodic boundary conditions";
}

// Test 2: 90-degree rotation
TEST_F(MolecularDatabaseTest, Rotation90Degrees) {
    mol2.apos[0].set(0.0, 0.0, 0.0);
    mol2.apos[1].set(0.0, -1.0, 0.0);
    mol2.apos[2].set(1.0, -1.0, 0.0);
    mol2.apos[3].set(1.0, 0.0, 0.0);
    
    double distance = db.computeDistanceOnSurf(&mol1, &mol2, &lat_vec);
    EXPECT_LT(distance, 0.1) << "Distance between 90-degree rotated molecules should be near zero";
}

// Test 3: 45-degree rotation
TEST_F(MolecularDatabaseTest, Rotation45Degrees) {
    mol2.apos[0].set(0.0, 0.0, 0.0);
    mol2.apos[1].set(0.7, 0.7, 0.0);
    mol2.apos[2].set(1.4, 0.0, 0.0);
    mol2.apos[3].set(0.7, -0.7, 0.0);
    
    double distance = db.computeDistanceOnSurf(&mol1, &mol2, &lat_vec);
    EXPECT_GT(distance, 0.1) << "Distance between 45-degree rotated molecules should be greater than threshold";
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
