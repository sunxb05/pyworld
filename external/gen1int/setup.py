#!/usr/bin/env python
"""Simple Python building script for Gen1Int."""

__author__ = "Bin Gao, and Andreas J. Thorvaldsen"
__copyright__ = "Copyright 2009-2012"
__credits__ = ["Radovan Bast", "Kenneth Ruud"]
__license__ = "LGPLv3"
__version__ = "0.2.1"
__maintainer__ = "Bin Gao"
__email__ = "bin.gao@uit.no"
__status__ = "Development"

gen1int_description = """\
Gen1Int is a Fortran 90 library (with Python interface) to evaluate the
derivatives of one-electron integrals with respect to the geometry perturbation,
external electric and magnetic fields, and total rotational angular momentum
at zero fields with contracted rotational London atomic orbitals (LAO)."""

gen1int_classifiers = [\
'Development Status :: 3 - Alpha',\
'Intended Audience :: by End-User Class :: Developers',\
'Intended Audience :: by Industry or Sector :: Science/Research',\
'License :: OSI-Approved Open Source :: GNU Library or "Lesser" General Public License version 3.0 (LGPLv3)',\
'Topic :: Scientific/Engineering :: Chemistry',\
'Topic :: Scientific/Engineering :: Molecular Science',\
'Topic :: Scientific/Engineering :: Physics',\
'Programming Language :: C',\
'Programming Language :: Fortran90',\
'Programming Language :: Python',\
'Operating System :: Microsoft :: Windows',\
'Operating System :: Unix',\
'Operating System :: MacOS',\
'User Interface :: Textual :: Console/Terminal']

# Preprocessed source codes of Gen1Int.Tools
fpp_tools = ['src/tools/py_norm_contr_cgto.F90', \
             'src/tools/py_norm_contr_sgto.F90', \
             'src/tools/py_reorder_ints.F90',    \
             'src/tools/py_trace_ints.F90']
# Preprocessed source codes of Gen1Int.ContrInt
fpp_contrint = ['src/py_error_stop.F90',                \
                'src/basic/py_binom_coeff.F90',         \
                'src/basic/py_const_contr_ints.F90',    \
                'src/basic/py_hgto_to_cgto.F90',        \
                'src/basic/py_hgto_to_sgto.F90',        \
                'src/basic/py_shell_scatter.F90',       \
                'src/basic/py_sort_cents.F90',          \
                'src/basic/py_next_permutation.F90',    \
                'src/auxfun/py_aux_boys_vec.F90',       \
                'src/geom/py_geom_total.F90',           \
                'src/geom/py_geom_part_zero.F90',       \
                'src/geom/py_geom_part_one.F90',        \
                'src/mag/py_hgto_to_lcgto.F90',         \
                'src/mag/py_london_mom_hgto.F90',       \
                'src/carmom/py_carmom_deriv.F90',       \
                'src/carmom/py_carmom_hbra.F90',        \
                'src/carmom/py_carmom_hrr_ket.F90',     \
                'src/carmom/py_carmom_moment.F90',      \
                'src/carmom/py_prim_hgto_carmom.F90',   \
                'src/carmom/py_contr_cgto_carmom.F90',  \
                'src/carmom/py_contr_sgto_carmom.F90',  \
                'src/carmom/py_contr_csgto_carmom.F90', \
                'src/delta/py_delta_geom.F90',          \
                'src/delta/py_delta_hket.F90',          \
                'src/delta/py_delta_moment.F90',        \
                'src/delta/py_prim_hgto_delta.F90',     \
                'src/delta/py_contr_cgto_delta.F90',    \
                'src/delta/py_contr_sgto_delta.F90',    \
                'src/nucpot/py_nucpot_geom.F90',        \
                'src/nucpot/py_nucpot_hket.F90',        \
                'src/nucpot/py_nucpot_hbra.F90',        \
                'src/nucpot/py_prim_hgto_nucpot.F90',   \
                'src/nucpot/py_contr_cgto_nucpot.F90',  \
                'src/nucpot/py_contr_sgto_nucpot.F90',  \
                'src/gaupot/py_gaupot_geom.F90',        \
                'src/gaupot/py_prim_hgto_gaupot.F90',   \
                'src/gaupot/py_contr_cgto_gaupot.F90',  \
                'src/gaupot/py_contr_sgto_gaupot.F90',  \
                'src/odist/py_prim_hgto_odist.F90',     \
                'src/odist/py_contr_cgto_odist.F90',    \
                'src/odist/py_contr_sgto_odist.F90',    \
                'src/value/py_prim_hgto_value.F90',     \
                'src/value/py_const_contr_gto.F90',     \
                'src/value/py_contr_cgto_value.F90',    \
                'src/value/py_contr_sgto_value.F90']
# Defined preprocessor directive options
def_opts = {'DEBUG':False,'XTIME':False}

# Gets the name of header file in the line of #include
def getHeader(src_line):
    # Case of #include "foobar.h"
    split_line = src_line.split('"')
    if len(split_line)>=2:
        inc_file = split_line[1]
    else:
        # Case of #include 'foobar.h'
        try:
            split_line = src_line.split("'")
            inc_file = split_line[1]
        # Case of #include <foobar.h>
        except:
            split_line = src_line.split('<')
            inc_file = split_line[1].split('>')[0]
        finally:
            raise Exception('Can not find include file in '+src_line)
    return inc_file

# Substitutes the parameters in #define
def subDefine(src_line, def_param):
    out_str = src_line
    # Replaced by predefined directives
    for def_key in def_param.keys():
        tmp_str = out_str.replace(def_key, def_param[def_key])
        out_str = tmp_str
    return out_str

# Takes care #ifdef preprocessor directive
def preIfDef(src_lines, def_param):
    key_opt = src_lines[0].split()[1]
    out_str = ''
    num_lines = len(src_lines)
    idx_line = 1
    # Dumps the following statements till #else
    if def_opts.get(key_opt):
        while idx_line<num_lines:
            if src_lines[idx_line][0:5] == '#else':
                idx_line += 1
                while idx_line<num_lines:
                    if src_lines[idx_line][0:6] == '#endif':
                        break
                    idx_line += 1
                break
            elif src_lines[idx_line][0:6] == '#endif':
                break
            # Another #ifdef preprocessor directive
            elif src_lines[idx_line][0:6] == '#ifdef':
                [tmp_str,num_pass_lines] = preIfDef(src_lines[idx_line:], def_param)
                out_str += tmp_str
                idx_line += num_pass_lines
            # Another #ifndef preprocessor directive
            elif src_lines[idx_line][0:7] == '#ifndef':
                [tmp_str,num_pass_lines] = preIfNDef(src_lines[idx_line:], def_param)
                out_str += tmp_str
                idx_line += num_pass_lines
            # Another #if preprocessor directive
            elif src_lines[idx_line][0:3] == '#if':
                [tmp_str,num_pass_lines] = preIfBlock(src_lines[idx_line:], def_param)
                out_str += tmp_str
                idx_line += num_pass_lines
            # Another #include preprocessor directive
            elif src_lines[idx_line][0:8] == '#include':
                inc_file = getHeader(src_lines[idx_line][8:])
                out_str += preSrc('include/'+inc_file, def_param)
                idx_line += 1
            # Another #define preprocessor directive
            elif src_lines[idx_line][0:7] == '#define':
                split_line = src_lines[idx_line].split()
                # By considering already defined parameters, such as REALK
                val_param = subDefine(split_line[2], def_param)
                def_param.update({split_line[1]:val_param})
                idx_line += 1
            else:
                out_str += subDefine(src_lines[idx_line], def_param)
                idx_line += 1
    # Dumps the statements after #else if there is
    else:
        while idx_line<num_lines:
            if src_lines[idx_line][0:5] == '#else':
                idx_line += 1
                while idx_line<num_lines:
                    if src_lines[idx_line][0:6] == '#endif':
                        break
                    # Another #ifdef preprocessor directive
                    elif src_lines[idx_line][0:6] == '#ifdef':
                        [tmp_str,num_pass_lines] = preIfDef(src_lines[idx_line:], def_param)
                        out_str += tmp_str
                        idx_line += num_pass_lines
                    # Another #ifndef preprocessor directive
                    elif src_lines[idx_line][0:7] == '#ifndef':
                        [tmp_str,num_pass_lines] = preIfNDef(src_lines[idx_line:], def_param)
                        out_str += tmp_str
                        idx_line += num_pass_lines
                    # Another #if preprocessor directive
                    elif src_lines[idx_line][0:3] == '#if':
                        [tmp_str,num_pass_lines] = preIfBlock(src_lines[idx_line:], def_param)
                        out_str += tmp_str
                        idx_line += num_pass_lines
                    # Another #include preprocessor directive
                    elif src_lines[idx_line][0:8] == '#include':
                        inc_file = getHeader(src_lines[idx_line][8:])
                        out_str += preSrc('include/'+inc_file, def_param)
                        idx_line += 1
                    # Another #define preprocessor directive
                    elif src_lines[idx_line][0:7] == '#define':
                        split_line = src_lines[idx_line].split()
                        # By considering already defined parameters, such as REALK
                        val_param = subDefine(split_line[2], def_param)
                        def_param.update({split_line[1]:val_param})
                        idx_line += 1
                    else:
                        out_str += subDefine(src_lines[idx_line], def_param)
                        idx_line += 1
                break
            elif src_lines[idx_line][0:6] == '#endif':
                break
            else:
                idx_line += 1
    return out_str,idx_line+1

# Takes care #ifndef preprocessor directive
def preIfNDef(src_lines, def_param):
    key_opt = src_lines[0].split()[1]
    out_str = ''
    num_lines = len(src_lines)
    idx_line = 1
    # Dumps the statements after #else if there is
    if def_opts.get(key_opt):
        while idx_line<num_lines:
            if src_lines[idx_line][0:5] == '#else':
                idx_line += 1
                while idx_line<num_lines:
                    if src_lines[idx_line][0:6] == '#endif':
                        break
                    # Another #ifdef preprocessor directive
                    elif src_lines[idx_line][0:6] == '#ifdef':
                        [tmp_str,num_pass_lines] = preIfDef(src_lines[idx_line:], def_param)
                        out_str += tmp_str
                        idx_line += num_pass_lines
                    # Another #ifndef preprocessor directive
                    elif src_lines[idx_line][0:7] == '#ifndef':
                        [tmp_str,num_pass_lines] = preIfNDef(src_lines[idx_line:], def_param)
                        out_str += tmp_str
                        idx_line += num_pass_lines
                    # Another #if preprocessor directive
                    elif src_lines[idx_line][0:3] == '#if':
                        [tmp_str,num_pass_lines] = preIfBlock(src_lines[idx_line:], def_param)
                        out_str += tmp_str
                        idx_line += num_pass_lines
                    # Another #include preprocessor directive
                    elif src_lines[idx_line][0:8] == '#include':
                        inc_file = getHeader(src_lines[idx_line][8:])
                        out_str += preSrc('include/'+inc_file, def_param)
                        idx_line += 1
                    # Another #define preprocessor directive
                    elif src_lines[idx_line][0:7] == '#define':
                        split_line = src_lines[idx_line].split()
                        # By considering already defined parameters, such as REALK
                        val_param = subDefine(split_line[2], def_param)
                        def_param.update({split_line[1]:val_param})
                        idx_line += 1
                    else:
                        out_str += subDefine(src_lines[idx_line], def_param)
                        idx_line += 1
                break
            elif src_lines[idx_line][0:6] == '#endif':
                break
            else:
                idx_line += 1
    # Dumps the following statements till #else
    else:
        while idx_line<num_lines:
            if src_lines[idx_line][0:5] == '#else':
                idx_line += 1
                while idx_line<num_lines:
                    if src_lines[idx_line][0:6] == '#endif':
                        break
                    idx_line += 1
                break
            elif src_lines[idx_line][0:6] == '#endif':
                break
            # Another #ifdef preprocessor directive
            elif src_lines[idx_line][0:6] == '#ifdef':
                [tmp_str,num_pass_lines] = preIfDef(src_lines[idx_line:], def_param)
                out_str += tmp_str
                idx_line += num_pass_lines
            # Another #ifndef preprocessor directive
            elif src_lines[idx_line][0:7] == '#ifndef':
                [tmp_str,num_pass_lines] = preIfNDef(src_lines[idx_line:], def_param)
                out_str += tmp_str
                idx_line += num_pass_lines
            # Another #if preprocessor directive
            elif src_lines[idx_line][0:3] == '#if':
                [tmp_str,num_pass_lines] = preIfBlock(src_lines[idx_line:], def_param)
                out_str += tmp_str
                idx_line += num_pass_lines
            # Another #include preprocessor directive
            elif src_lines[idx_line][0:8] == '#include':
                inc_file = getHeader(src_lines[idx_line][8:])
                out_str += preSrc('include/'+inc_file, def_param)
                idx_line += 1
            # Another #define preprocessor directive
            elif src_lines[idx_line][0:7] == '#define':
                split_line = src_lines[idx_line].split()
                # By considering already defined parameters, such as REALK
                val_param = subDefine(split_line[2], def_param)
                def_param.update({split_line[1]:val_param})
                idx_line += 1
            else:
                out_str += subDefine(src_lines[idx_line], def_param)
                idx_line += 1
    return out_str,idx_line+1

# Takes care #if block
def preIfBlock(src_lines, def_param):
    # Processes the #if statement
    if_cond = src_lines[0][4:].replace('\n','')
    for def_key in def_param.keys():
        tmp_str = if_cond.replace(def_key, def_param[def_key])
        if_cond = tmp_str
    for def_key in def_opts.keys():
        tmp_str = if_cond.replace(def_key, repr(def_opts[def_key]))
        if_cond = tmp_str
    tmp_str = if_cond.replace('defined','')
    if_cond = tmp_str.replace(' ','')
    tmp_str = if_cond.replace('!(True)','False')
    if_cond = tmp_str.replace('!(False)','True')
    # Maximum number of iterations, for safety
    max_num_iter = 900
    # Starts to decode the #if statement
    iter = 0
    while if_cond!='True' and if_cond!='False' and iter<max_num_iter:
        iter += 1
        tmp_str = if_cond.replace('(True)','True')
        if_cond = tmp_str.replace('(False)','False')
        tmp_str = if_cond.replace('True&&True','True')
        if_cond = tmp_str.replace('False&&True','False')
        tmp_str = if_cond.replace('True&&False','False')
        if_cond = tmp_str.replace('False&&False','False')
        tmp_str = if_cond.replace('True||True','True')
        if_cond = tmp_str.replace('False||True','True')
        tmp_str = if_cond.replace('True||False','True')
        if_cond = tmp_str.replace('False||False','False')
    # Checks the result
    if if_cond=='True':
        dump_if_block = True
    elif if_cond=='False':
        dump_if_block = False
    else:
        raise StopIteration('Can not process #if statement '+src_lines[0][4:])
    # Processes the #if block
    out_str = ''
    num_lines = len(src_lines)
    idx_line = 1
    # Dumps the following statements till #else
    if dump_if_block:
        while idx_line<num_lines:
            if src_lines[idx_line][0:5] == '#else':
                idx_line += 1
                while idx_line<num_lines:
                    if src_lines[idx_line][0:6] == '#endif':
                        break
                    idx_line += 1
                break
            elif src_lines[idx_line][0:6] == '#endif':
                break
            # Another #ifdef preprocessor directive
            elif src_lines[idx_line][0:6] == '#ifdef':
                [tmp_str,num_pass_lines] = preIfDef(src_lines[idx_line:], def_param)
                out_str += tmp_str
                idx_line += num_pass_lines
            # Another #ifndef preprocessor directive
            elif src_lines[idx_line][0:7] == '#ifndef':
                [tmp_str,num_pass_lines] = preIfNDef(src_lines[idx_line:], def_param)
                out_str += tmp_str
                idx_line += num_pass_lines
            # Another #if preprocessor directive
            elif src_lines[idx_line][0:3] == '#if':
                [tmp_str,num_pass_lines] = preIfBlock(src_lines[idx_line:], def_param)
                out_str += tmp_str
                idx_line += num_pass_lines
            # Another #include preprocessor directive
            elif src_lines[idx_line][0:8] == '#include':
                inc_file = getHeader(src_lines[idx_line][8:])
                out_str += preSrc('include/'+inc_file, def_param)
                idx_line += 1
            # Another #define preprocessor directive
            elif src_lines[idx_line][0:7] == '#define':
                split_line = src_lines[idx_line].split()
                # By considering already defined parameters, such as REALK
                val_param = subDefine(split_line[2], def_param)
                def_param.update({split_line[1]:val_param})
                idx_line += 1
            else:
                out_str += subDefine(src_lines[idx_line], def_param)
                idx_line += 1
    # Dumps the statements after #else if there is
    else:
        while idx_line<num_lines:
            if src_lines[idx_line][0:5] == '#else':
                idx_line += 1
                while idx_line<num_lines:
                    if src_lines[idx_line][0:6] == '#endif':
                        break
                    # Another #ifdef preprocessor directive
                    elif src_lines[idx_line][0:6] == '#ifdef':
                        [tmp_str,num_pass_lines] = preIfDef(src_lines[idx_line:], def_param)
                        out_str += tmp_str
                        idx_line += num_pass_lines
                    # Another #ifndef preprocessor directive
                    elif src_lines[idx_line][0:7] == '#ifndef':
                        [tmp_str,num_pass_lines] = preIfNDef(src_lines[idx_line:], def_param)
                        out_str += tmp_str
                        idx_line += num_pass_lines
                    # Another #if preprocessor directive
                    elif src_lines[idx_line][0:3] == '#if':
                        [tmp_str,num_pass_lines] = preIfBlock(src_lines[idx_line:], def_param)
                        out_str += tmp_str
                        idx_line += num_pass_lines
                    # Another #include preprocessor directive
                    elif src_lines[idx_line][0:8] == '#include':
                        inc_file = getHeader(src_lines[idx_line][8:])
                        out_str += preSrc('include/'+inc_file, def_param)
                        idx_line += 1
                    # Another #define preprocessor directive
                    elif src_lines[idx_line][0:7] == '#define':
                        split_line = src_lines[idx_line].split()
                        # By considering already defined parameters, such as REALK
                        val_param = subDefine(split_line[2], def_param)
                        def_param.update({split_line[1]:val_param})
                        idx_line += 1
                    else:
                        out_str += subDefine(src_lines[idx_line], def_param)
                        idx_line += 1
                break
            elif src_lines[idx_line][0:6] == '#endif':
                break
            else:
                idx_line += 1
    return out_str,idx_line+1

# Preprocesses Fortran 90 source code
def preSrc(src_code, def_param):
    fl_src = file(src_code)
    src_lines = fl_src.readlines()
    out_str = ''
    num_lines = len(src_lines)
    idx_line = 0
    while idx_line<num_lines:
        # Preprocessor directive #include
        if src_lines[idx_line][0:8] == '#include':
            inc_file = getHeader(src_lines[idx_line][8:])
            out_str += preSrc('include/'+inc_file, def_param)
            idx_line += 1
        # Preprocessor directive #ifdef
        elif src_lines[idx_line][0:6] == '#ifdef':
            [tmp_str,num_pass_lines] = preIfDef(src_lines[idx_line:], def_param)
            out_str += tmp_str
            idx_line += num_pass_lines
        # Preprocessor directive #ifndef
        elif src_lines[idx_line][0:7] == '#ifndef':
            [tmp_str,num_pass_lines] = preIfNDef(src_lines[idx_line:], def_param)
            out_str += tmp_str
            idx_line += num_pass_lines
        # Preprocessor directive #if
        elif src_lines[idx_line][0:3] == '#if':
            [tmp_str,num_pass_lines] = preIfBlock(src_lines[idx_line:], def_param)
            out_str += tmp_str
            idx_line += num_pass_lines
        # Preprocessor directive #define
        elif src_lines[idx_line][0:7] == '#define':
            split_line = src_lines[idx_line].split()
            # By considering already defined parameters, such as REALK
            val_param = subDefine(split_line[2], def_param)
            def_param.update({split_line[1]:val_param})
            idx_line += 1
        else:
            out_str += subDefine(src_lines[idx_line], def_param)
            idx_line += 1
    fl_src.close()
    return out_str

# Removes files
def rmFiles(list_files=[]):
    from os import remove
    for ifile in list_files:
        remove(ifile)

# Configuration constructor
def setConfig(parent_package="", top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(None, parent_package, top_path)
    config.add_subpackage("Gen1Int")
    config.add_extension('Gen1Int.Tools', sources=fpp_tools, \
                         f2py_options=['--no-wrap-functions'])
    config.add_extension('Gen1Int.ContrInt', sources=fpp_contrint, \
                         f2py_options=['--no-wrap-functions'])
    return config

# Setup Gen1Int
def setupPackage():
    from numpy.distutils.core import setup
    # Preprocesses source codes
    for isrc in xrange(len(fpp_tools)):
        # Source code
        src_tools = "".join(fpp_tools[isrc].split('py_'))
        # Defined preprocessor directive parameters
        def_param = {}
        try:
            fl_fpp = file(fpp_tools[isrc], 'w')
            fl_fpp.write(preSrc(src_tools, def_param))
            fl_fpp.close()
        except:
            rmFiles(list_files=fpp_tools[0:isrc+1])
            raise Exception("Fatal error when preprocessing "+src_tools+"!\n")
    for isrc in xrange(len(fpp_contrint)):
        # Source code
        src_contrint = "".join(fpp_contrint[isrc].split('py_'))
        # Defined preprocessor directive parameters
        def_param = {}
        try:
            fl_fpp = file(fpp_contrint[isrc], 'w')
            fl_fpp.write(preSrc(src_contrint, def_param))
            fl_fpp.close()
        except:
            rmFiles(list_files=fpp_tools)
            rmFiles(list_files=fpp_contrint[0:isrc+1])
            raise Exception("Fatal error when preprocessing "+src_contrint+"!\n")
    # Setup Gen1Int
    setup(
        name="Gen1Int",
        version="0.2.1",
        author="Bin Gao, and Andreas J. Thorvaldsen",
        author_email="bin.gao@uit.no, and andreas.thorvaldsen@uit.no",
        maintainer="Bin Gao",
        maintainer_email="bin.gao@uit.no",
        url="http://repo.ctcc.no/projects/gen1int",
        download_url="http://sourceforge.net/projects/gen1int",
        license="LGPLv3",
        description="Generalized One-Electron Integral library",
        long_description=gen1int_description,
        platforms=['Windows','Unix','MacOS'],
        classifiers=gen1int_classifiers,
        configuration=setConfig)
    # Removes the preprocessing source codes
    rmFiles(list_files=fpp_tools)
    rmFiles(list_files=fpp_contrint)
    return

if __name__ == '__main__':
    setupPackage()

#ScipyTest('Gen1Int').test(level=1,verbosity=1)
