# Copyright (C) 2013 - 2019 by the authors of the ASPECT code.
#
# This file is part of ASPECT.
#
# ASPECT is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# ASPECT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ASPECT; see the file doc/COPYING.  If not see
# <http://www.gnu.org/licenses/>.

project(opendap-sph Fortran)

add_executable(opendap_convert opendap_convert.f)
add_executable(sphadd sphadd.f)

target_include_directories(opendap_convert PUBLIC
        /usr/include)

target_include_directories(sphadd PUBLIC
        /usr/include)
