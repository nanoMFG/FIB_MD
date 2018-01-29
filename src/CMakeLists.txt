#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.4)

project(FIB_MD-src Fortran C)

set ( srcs 
  ${CMAKE_CURRENT_SOURCE_DIR}/atom_traj.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/bubble_dim_life.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/count_sputter.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/data.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/hole_den_time.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/hole_profile2.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/hole_profile3.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/hole_profile.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/ion_traj.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/ion_traj_from_MD.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/knock.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/liquidprof.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/main.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/meanvel_relax.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/midsection.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/nlist.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/out.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/outfiles.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/parallel.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/par.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/prms.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/ran1.f
  ${CMAKE_CURRENT_SOURCE_DIR}/rhorms.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/si.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/solidprof.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/stats.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/temp.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/test.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/testprog.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/timeint.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/tools.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/tsi.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/upsidedown_den_of_crystal.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/upsidedown_disp.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/upsidedown.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/upsidedown_para.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/velocity_avg.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/void_count.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/write_sput_files2.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/write_sput_files.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/write_XZ_plane.f90
)

add_executable (mdrun2 ${srcs})
