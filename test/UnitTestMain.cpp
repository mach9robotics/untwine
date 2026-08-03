/*****************************************************************************
 *   Copyright (c) 2026, Mach9 Robotics                                       *
 *                                                                           *
 *   Entry point shared by all of untwine's pure-logic unit + integration     *
 *   tests (everything except the end-to-end binary test in Tests.cpp).      *
 ****************************************************************************/

#include <gtest/gtest.h>

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
