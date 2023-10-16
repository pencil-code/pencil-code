from asyncore import write
import math
import re
import os
import csv
import argparse
import glob


def split_by_indexes(line,indexes):
    return [line[i:j] for i,j in zip(indexes, indexes[1:]+[None])]

def split_by_ops(line,op_chars):
    num_of_right_brackets = 0
    num_of_left_brackets = 0
    possible_var = ""
    split_indexes = [0]
    in_array = False
    for char_index,char in enumerate(line):
        if char == "(":
            num_of_left_brackets += 1
            if line[char_index-1]:
                back_index = char_index-1
                while(line[back_index] == " "):
                    back_index -= 1
                if line[back_index] not in "*-+-/(":
                    in_array = True
            else:
                if line[char_index-1] not in "*-+-/(":
                    # in_array = possible_var in local_variables or possible_var in self.static_variables 
                    in_array = True
            possible_var = ""
        elif char == ")":
            num_of_right_brackets += 1
            if num_of_left_brackets == num_of_right_brackets:
                in_array = False
            possible_var = ""
        elif (char in op_chars) and not in_array:
            split_indexes.append(char_index)
            possible_var = ""
        else:
            possible_var = possible_var + char
    res = [x for x in split_by_indexes(line,split_indexes) if x != ""]
    for i,val in enumerate(res):
        if val[0] in "+-*/":
            res[i] = val[1:]
    return [x for x in res if x != ""]
    
def okay_stencil_index(index,i):
    if index == ":":
        return True
    if i == 0:
        return index in ["l1:l2"]
    if i==1:
        return index in ["m"]
    if i==2:
        return index in ["n"]
def is_vector_stencil_index(index):
    if ":" not in index:
        return False
    lower,upper= [part.strip() for part in index.split(":")]
    if "(" in lower:
        lower_indexes = get_indexes(lower,lower.split("(",1)[0].strip(),0)
        upper_indexes= get_indexes(lower,lower.split("(",1)[0].strip(),0)
        lower = lower.split("(",1)[0].strip()
        upper = upper.split("(",1)[0].strip()
        assert(lower_indexes == upper_indexes)
        assert(len(lower_indexes) == 1)
        return lower[-1] == "x" and upper[-1] == "z"

    return lower[-1] == "x" and upper[-1] == "z"
def get_vtxbuf_name_from_index(prefix, index):
    ## $ will be replaced later when producing vector lines
    if ":" in index:
        lower,upper= [part.strip() for part in index.split(":")]
        if "(" in lower:
            lower_indexes = get_indexes(lower,lower.split("(",1)[0].strip(),0)
            upper_indexes= get_indexes(lower,lower.split("(",1)[0].strip(),0)
            assert(lower_indexes == upper_indexes)
            assert(len(lower_indexes) == 1)
            lower = lower.split("(",1)[0].strip()
            upper = upper.split("(",1)[0].strip()
            return f"{prefix}{(lower[1:-1]+lower_indexes[0]).upper()}$"
        return f"{prefix}{lower[1:-1].upper()}$"
    if "(" in index:
        indexes = get_indexes(index, index.split("(",1)[0].strip(),0)
        assert(len(indexes) == 1)
        index = index.split("(",1)[0].strip()
        return f"{prefix}{(index[1:]+indexes[0]).upper()}$"
    return f"{prefix}{index[1:].upper()}"

def get_function_call_index(function_call, lines):
    for i, (line,count) in enumerate(lines):
        if f"{function_call}(" in line or ("call" in line and function_call in line):
            return i
    print(f"Didn't find index for function_call: {function_call}")
    exit()
def has_balanced_parens(line):
    return sum([char == "(" for char in line]) == sum([char == ")" for char in line])
def replace_exp_once(line):
    num_of_left_brackets = 0
    num_of_right_brackets = 0
    starting_index = 0
    it_index = 1
    while starting_index == 0:
        if line[it_index] == "*" and line[it_index-1] == "*":
            starting_index = it_index-1
        it_index += 1
    forward_index = starting_index + 2
    backward_index = starting_index - 1
    save_index = starting_index - 1
    exponent = ""
    num_of_left_brackets == 0
    num_of_right_brackets == 0
    parse = True
    while parse:
        if line[forward_index] == "(":
            num_of_left_brackets += 1
        if line[forward_index] == ")":
            num_of_right_brackets += 1
        if (forward_index == len(line)-1 or line[forward_index+1] in " *+-;[]()/\{}") and num_of_left_brackets == num_of_right_brackets:
            parse = False
            exponent = exponent + line[forward_index]
        if parse:
            exponent = exponent + line[forward_index]
            forward_index += 1
    base = ""
    num_of_left_brackets == 0
    num_of_right_brackets == 0
    parse = True
    while parse:
        if line[backward_index] == "(":
            num_of_left_brackets += 1
        if line[backward_index] == ")":
            num_of_right_brackets += 1
        if (backward_index == 0 or line[backward_index-1] in " *+-;[]()/\{}")  and num_of_left_brackets == num_of_right_brackets:
            parse = False
            base = line[backward_index] + base
        if parse:
            base = line[backward_index] + base
            backward_index -= 1
    res = ""
    for i,x in enumerate(line):
        if i<backward_index or i>forward_index:
            res = res + x
        elif i ==backward_index:
            res = res + f"pow({base},{exponent})"
    return res
def replace_exp(line):
    while "**" in line:
        line = replace_exp_once(line)
    return line
def translate_to_c(type):
    if type =="real":
        return "AcReal"
    elif type =="logical":
        return "bool"
    elif type =="integer":
        return "int"
    elif type=="character":
        return "char*"
    else:
        print("WHAT is the translation for",type,"?")
        exit()
def common_data(list1, list2):
    result = False
 
    # traverse in the 1st list
    for x in list1:
 
        # traverse in the 2nd list
        for y in list2:
   
            # if one common
            if x == y:
                result = True
                return result
                 
    return result
def get_segment_indexes(segment,line,dims):
    return get_indexes(line[segment[1]:segment[2]],segment[0],dims)
def get_indexes(segment,var,dim):
    index_search = re.search(f"{var}\((.*)\)",segment)
    indexes = [":" for i in range(dim)]
    if index_search:
        indexes = []
        index = ""
        index_line = index_search.group(1)
        num_of_left_brackets = 0
        num_of_right_brackets = 0
        acc_str =""
        for char in index_line:
            acc_str = acc_str + char
            if char == "(": 
                num_of_left_brackets = num_of_left_brackets + 1
            if char == ")":
                num_of_left_brackets = num_of_right_brackets + 1
            if char == "," and num_of_left_brackets == num_of_right_brackets:
                indexes.append(index.split("(")[-1].split(")")[0].strip())
                index = ""
            else:
                index = index + char
        if index.strip() != "" and sum([char == "(" for char in index.strip()]) == sum([char == ")" for char in index.strip()]):
            indexes.append(index.strip())
        else:
            index = index[:-1]
            if index.strip() != "" and sum([char == "(" for char in index.strip()]) == sum([char == ")" for char in index.strip()]):
                indexes.append(index.strip())
            elif index.strip() != "":
                indexes.append(index.split("(")[-1].split(")")[0].strip())
    return indexes
def is_body_line(line):
    return not "::" in line and "subroutine" not in line and not line.split(" ")[0].strip() == "use" and "function" not in line
def is_init_line(line):
    return "::" in line
def add_splits(line):
    res = ""
    for line_part in line.split("\n"):
        res = res + add_splits_per_line(line_part) + "\n"
    return res
def add_splits_per_line(line):
    res = ""
    iterator = 0
    split_modulus = 80
    for char in line:
        res = res + char
        iterator = iterator + 1
        if char == " " and iterator >=split_modulus:
            iterator = 0
            res = res + "&\n"
    return res



def is_use_line(line):
    return "use" == line.split(" ")[0].strip()
def is_declaration_line(line):
    return "::" in line
def new_is_variable_line(line):
    parts = line.split("::")
    if "!" in parts[0]:
        return False
    return len(parts)>1
def is_variable_line(line_elem):
    parts = line_elem[0].split("::")
    if "!" in parts[0]:
        return False
    return len(parts)>1
def merge_dictionaries(dict1, dict2):
    merged_dict = dict1.copy()
    merged_dict.update(dict2)
    return merged_dict
def get_variable_segments(line, variables):
    check_string = ""
    accumulating = False
    start_index = 0
    var = ""
    res = []
    i = 0
    inside_indexes = False
    num_of_left_brackets = 0
    num_of_right_brackets = 0
    while i < len(line):
        char = line[i]
        if inside_indexes and char == "(":
            num_of_left_brackets += 1
        if inside_indexes and char == ")":
            num_of_right_brackets += 1
            if  num_of_left_brackets == num_of_right_brackets:
                inside_indexes = False
        if char in " /*+-<>=()" and not inside_indexes:
            if accumulating:
                accumulating = False
                if char == ")":
                    res.append((var,start_index,i+1))
                    i = i +1
                else:
                    res.append((var,start_index,i))
            check_string = ""
            start_index = i + 1
        else:
            check_string = check_string + char
        if check_string in variables:
            if i >= len(line)-1:
                res.append((check_string,start_index,i+1))
                return res
            elif line[i+1] in " /*+-<>=()":
                var = check_string
                if line[i+1] != "(":
                    res.append((var,start_index,i+1))
                else:
                    inside_indexes = True
                    num_of_left_brackets = 1
                    num_of_right_brackets = 0
                    accumulating = True
                    i = i + 1
        i = i + 1
    return res
def replace_variable(line, old_var, new_var):
    check_string = ""
    indexes = []
    for i, char in enumerate(line):
        if char in  " ,/*+-<>=()[];:.":
            check_string = ""
        else:
            check_string = check_string + char
        if check_string == old_var:
            if i == len(line)-1:
                indexes.append((i+1-len(old_var),i+1))
            elif line[i+1] in   " ,/*+-<>=()[];:.":
                indexes.append((i+1-len(old_var),i+1))
    if len(indexes) == 0:
        return line
    res_line = ""
    last_index = 0
    for lower,upper in indexes:
        res_line= res_line+ line[last_index:lower]
        res_line= res_line+ new_var
        last_index = upper
    res_line = res_line + line[indexes[-1][1]:]
    search_var = new_var.split("(")[0].strip()
    new_var_segs = get_variable_segments(res_line,[search_var])
    if len(new_var_segs) == 0:
        return res_line
    seg_out = [res_line[seg[1]:seg[2]] for seg in new_var_segs]
    seg_indexes = [(seg[1],seg[2]) for seg in new_var_segs]
    for seg_index, seg in enumerate(new_var_segs):
        if seg[2] < len(res_line)-1:
            my_bool = res_line[seg[2]] == "(" and res_line[seg[2]-1] == ")"
            if my_bool:
                old_indexes = get_indexes(f"{res_line[seg[1]:seg[2]]})",search_var,-1)
                start_index = seg[2]+1
                end_index = start_index
                num_of_left_brackets = 0
                num_of_right_brackets = 0
                while res_line[end_index] != ")" or (num_of_left_brackets >= num_of_right_brackets + 1):
                    if res_line[end_index] == ")":
                        num_of_right_brackets = num_of_right_brackets + 1
                    if res_line[end_index] == "(":
                        num_of_left_brackets = num_of_left_brackets + 1
                    end_index = end_index + 1
                end_index = end_index + 1
                new_indexes = get_indexes(f"{search_var}({res_line[start_index:end_index]})",search_var,-1)
                combined_indexes = []
                new_indexes_iterator = 0 
                for old_index in old_indexes:
                    if ":" not in old_index:
                        combined_indexes.append(old_index)
                    elif old_index == ":":
                        combined_indexes.append(new_indexes[new_indexes_iterator])
                        new_indexes_iterator = new_indexes_iterator + 1
                    elif new_indexes[new_indexes_iterator] == ":":
                        combined_indexes.append(old_index)
                        new_indexes_iterator = new_indexes_iterator + 1
                    else:
                        old_lower, old_upper = [part.strip() for part in old_index.split(":")]
                        if ":" in new_indexes[new_indexes_iterator]:
                            new_lower, old_lower = [part.strip() for part in new_indexes[new_indexes_iterator].split(":")]
                            combined_indexes.append(f"{old_lower}+({new_lower}-1):{old_upper}-({new_upper}-1)")
                        else:
                            combined_indexes.append(f"{old_lower}+({new_indexes[new_indexes_iterator].strip()}-1)")
                        new_indexes_iterator = new_indexes_iterator + 1
                combined_indexes = [index.strip() for index in combined_indexes]
                sg_res = search_var + "(" 
                for i, index in enumerate(combined_indexes):
                    sg_res = sg_res + index
                    if i<len(combined_indexes)-1:
                        sg_res = sg_res + ","
                sg_res = sg_res + ")"
                seg_out[seg_index] = sg_res
                seg_indexes[seg_index] = (seg[1],end_index)
    last_res_line = ""
    last_index = 0
    for i,index in enumerate(seg_indexes):
        lower,upper = index
        last_res_line= last_res_line+ res_line[last_index:lower]
        last_res_line= last_res_line+ seg_out[i] 
        last_index = upper
    last_res_line = last_res_line + res_line[seg_indexes[-1][1]:]
    return last_res_line 

def second_val(x):
    return x[1]
def parse_variable(line_segment):
    iter_index = len(line_segment)-1
    end_index = iter_index
    start_index = 0
    num_of_left_brackets = 0
    num_of_right_brackets = 0
    while iter_index>0 and (line_segment[iter_index] in " ()" or num_of_left_brackets != num_of_right_brackets):
        elem = line_segment[iter_index]
        if elem == "(":
            num_of_left_brackets += 1
        elif elem == ")":
            num_of_right_brackets += 1
        iter_index -= 1
    end_index = iter_index+1

    while iter_index>0 and line_segment[iter_index] not in " *+-();":
        iter_index -= 1
        
    if iter_index == 0:
        res = line_segment[0:end_index]
    else:
        res = line_segment[iter_index+1:end_index]
    if "%" in res:
        end_index = 0
        while res[end_index] != "%":
            end_index += 1
        res = res[0:end_index]
    return res

def get_used_variables_from_line(line):
    characters_to_space= ["/",";",":", ",","+","-","*","(",")","=","<","!",">","[","]"]
    for character in characters_to_space:
        line = line.replace(character, " ")
    parts = [part.strip() for part in line.split(" ")]
    return [parse_variable(part) for part in parts if part]
def get_used_variables(lines):
    res = []
    for line in lines:
        res.extend(get_used_variables_from_line(line))
    return res
def get_var_segments_in_line(line,variables):
    vars = [var for var in get_used_variables_from_line(line) if var in variables]
    res = get_variable_segments(line, vars)
    return sorted(res, key = second_val)

def get_rhs_segment(line):
    res = []
    index = 0
    start_index = 0
    num_of_single_quotes = 0
    num_of_double_quotes = 0
    num_of_left_brackets = 0
    num_of_right_brackets = 0
    while index<len(line)-1:
        index += 1
        elem = line[index]
        if elem == "'":
            num_of_single_quotes += 1
        elif elem == '"':
            num_of_double_quotes += 1
        elif elem == "=" and line[index-1] not in "<>=!" and line[index+1] not in "<>=!" and num_of_single_quotes%2==0 and num_of_double_quotes%2==0 and num_of_left_brackets == num_of_right_brackets:
            return line[start_index:index]
def get_dims_from_indexes(indexes,rhs_var):
    dims = []
    num_of_looped_dims = 0
    for i,index in enumerate(indexes):
        if index == ":":
            dims.append(f"size({rhs_var},dim={i+1})")
            num_of_looped_dims = num_of_looped_dims + 1
        elif ":" in index:
            lower, upper = [var.strip() for var in index.split(":")]
            dims.append(f"{upper}-{lower}+1")
            num_of_looped_dims = num_of_looped_dims + 1
    return (dims, num_of_looped_dims)


def get_new_indexes(segment,var,dim):
    indexes = get_indexes(segment,var,dim) 
    transformed_indexes = 0
    new_indexes = []
    for index in indexes: 
        if index == ":":
            new_indexes.append(f"explicit_index_{transformed_indexes}")
            transformed_indexes = transformed_indexes + 1
        elif ":" in index:
            lower, upper = [var.strip() for var in index.split(":")]
            new_indexes.append(f"{lower} + explicit_index_{transformed_indexes}")
            transformed_indexes = transformed_indexes + 1
        else:
            new_indexes.append(index)
    return new_indexes
def build_new_access(var,new_indexes):
    res = f"{var}("
    for i,index in enumerate(new_indexes):
        res = res + index
        if i < len(new_indexes)-1:
            res = res + ","
    res = res + ")"
    return res


class Parser:

    def __init__(self, files,config):
        self.known_dims = {"beta_glnrho_scaled": ["3"]}
        self.test_to_c = config["to_c"]
        self.known_values = {}
        self.inline_num = 0
        self.sample_dir = config["sample_dir"]
        self.offload_type = None
        self.offloading = config["offload"]
        #we don't want print out values
        self.ranges = []
        self.flag_mappings= {
            "headtt":".false.",
            "ldebug":".false.",
        }
        self.default_mappings = {
            "ieos_profile": "'nothing'",
            "mass_source_profile":"'nothing'",
            "lschur_3d3d1d": ".false.",
            "lconservative": ".false.",
            "ldiff_hyper3_polar": ".true.",
            "ldiff_hyper3":".true.",
            "ldiff_hyper3_strict":".true.",
            "ldiff_hyper3_mesh": ".false.",
            "ldiff_hyper3_aniso" :".false.",
            "ldiff_hyper3lnrho": ".false.",
            "ldiff_hyper3lnrho_strict": ".false.",
            "lcylindrical_coords": ".false.",
            "lspherical_coords": ".false.",
            #The following looked up from source code
            "lcontinuity_gas": "(lhydro .or. lhydro_kinematic) .and. nwgrid/=1",
            "lupdate_mass_source": "lmass_source .and. t>=tstart_mass_source .and. (tstop_mass_source==-1.0 .or. t<=tstop_mass_source)",
            "llorentz_rhoref": "llorentzforce .and. rhoref/=impossible .and. rhoref>0."

        }
        self.safe_subs_to_remove = ["print","not_implemented","fatal_error","keep_compiler_quiet"]
        # self.profile_mappings = {
        #     "sinth": "PROFILE_Y", 
        #     "cosph": "PROFILE_Z",
        #     "costh":"PROFILE_Y",
        #     "sinph":"PROFILE_Z",
        #     "x":"PROFILE_X",
        #     "y":"PROFILE_Y",
        #     "z": "PROFILE_Z",
        #     "dline_1_x": "PROFILE_X",
        #     "dline_1_y": "PROFILE_Y",
        #     "dline_1_z": "PROFILE_Z",
        #     "etatotal": "PROFILE_X",
        #     "cs2cool_x": "PROFILE_X",
        #     "cs2mxy": "PROFILE_XY",
        #     "cs2mx": "PROFILE_X",
        #     }
        self.func_info = {}
        self.file_info = {}
        self.module_info = {}
        self.static_variables = {}
        self.lines = {}
        self.parsed_files_for_static_variables = []
        self.parsed_modules = []
        self.parsed_subroutines = []
        self.loaded_files = []
        self.subroutine_order = 0
        self.used_files = []
        self.found_function_calls = []
        self.used_files = files 
        self.used_static_variables = []
        self.functions_in_file = {}
        self.static_writes = []
        self.rewritten_functions = {}
        self.directory = config["directory"]
        self.subroutine_modifies_param = {}
        self.struct_table = {}
        ignored_files = []
        self.get_chosen_modules(config["Makefile"])
        self.ignored_modules = ["hdf5"]
        self.ignored_files = ["nodebug.f90","/boundcond_examples/","deriv_alt.f90","boundcond_alt.f90", "diagnostics_outlog.f90","pscalar.f90", "/cuda/", "/obsolete/", "/inactive/", "/astaroth/", "/initial_condition/", "/pre_and_post_processing/", "/scripts/"]
        self.ignored_files.append("magnetic_ffreeMHDrel.f90")
        self.ignored_files.append("photoelectric_dust.f90")
        self.ignored_files.append("interstellar_old.f90")
        self.ignored_files.append("spiegel.f90")
        self.used_files = [file for file in self.used_files if not any([ignored_file in file for ignored_file in self.ignored_files])  and ".f90" in file]
        for file in self.used_files:
            self.get_lines(file)
        self.used_files = [file for file in self.used_files if not self.file_info[file]["is_program_file"]]
        self.used_files.append(f"{self.directory}/run.f90")
        # self.used_files = list(filter(lambda x: x not in ignored_files and not self.is_program_file(x) and "/inactive/" not in x and "deriv_alt.f90" not in x and "diagnostics_outlog.f90" not in x "boundcond_alt.f90" and not re.search("cuda",x) and re.search("\.f90", x) and not re.search("astaroth",x)  and "/obsolete/" not in x,self.used_files))
        self.ignored_subroutines = ["alog10","count", "min1", "erf","aimag", "cmplx","len", "inquire", "floor", "matmul","ceiling", "achar", "adjustl", "index", "iabs","tiny","dble","float","nullify","associated","nint","open","close","epsilon","random_seed","modulo","nearest","xor","ishft","iand","ieor","ior","random_number","all","any","deallocate","cshift","allocated","allocate","case","real","int","complex","character","if","elseif","where","while","elsewhere","forall","maxval", "minval", "dot_product", "abs", "alog", "mod", "size",  "sqrt", "sum","isnan", "exp", "spread", "present", "trim", "sign","min","max","sin","cos","log","log10","tan","tanh","cosh","sinh","asin","acos","atan","atan2","write","read","char"]
        ##Ask Matthias about these
        self.ignored_subroutines.extend(["DCONST","fatal_error", "terminal_highlight_fatal_error","warning","caller","caller2", "coeffsx","coeffsy","coeffsz","r1i","sth1i","yh_","not_implemented","die","deri_3d","u_dot_grad_mat"])
        self.module_variables = {}
        self.ignored_subroutines.extend(['keep_compiler_quiet','keep_compiler_quiet_r', 'keep_compiler_quiet_r1d', 'keep_compiler_quiet_r2d', 'keep_compiler_quiet_r3d', 'keep_compiler_quiet_r4d', 'keep_compiler_quiet_p', 'keep_compiler_quiet_bc', 'keep_compiler_quiet_sl', 'keep_compiler_quiet_i', 'keep_compiler_quiet_i1d', 'keep_compiler_quiet_i81d', 'keep_compiler_quiet_i2d', 'keep_compiler_quiet_i3d', 'keep_compiler_quiet_l', 'keep_compiler_quiet_l1d', 'keep_compiler_quiet_c', 'keep_compiler_quiet_c1d',"calc_del6_for_upwind"])
        self.ignored_subroutines.extend(["helflux", "curflux_ds"])
        self.ignored_subroutines.extend(["cffti","cffti1","cfftf","cfftb","kx_fft","ky_fft"])
        self.ignored_subroutines.extend(["bc_aa_pot2"])
        self.ignored_subroutines.extend(["caller", "caller0", "caller1", "caller2", "caller3", "caller4", "caller5", "caller5_str5"])
        self.ignored_subroutines.append("output_penciled_scal_c")
        self.ignored_subroutines.append("output_penciled_vect_c")
        self.ignored_subroutines.append("output_pencil_scal")
        self.ignored_subroutines.append("output_pencil_vect")
        self.ignored_subroutines.append("loptest")
        self.ignored_subroutines.append("result")
        #omp subs
        self.ignored_subroutines.extend(["omp_get_thread_num"])
        self.modules_in_scope = {}
        self.file = config["file"]
        # for file in self.used_files:
        #     self.get_lines(file)
    def map_to_new_index(self,index,i,local_variables,line,possible_values=[],make_sure_indexes_safe=False):
        if ":" in index:
            self.ranges.append((index,i,line))
        if index == ":":
            if i==0:
                return "idx.x"
            elif i==1:
                return "idx.y"
            elif i==2:
                return "idx.z"
            else:
                print("whole range in last f index?")
                print(line)
                exit()
        lower, upper = [part.strip() for part in index.split(":")]
        if index == "idx.x":
            return "idx.x"
        if index == "idx.y":
            return "idx.y"
        if index == "idx.z":
            return "idx.z"

        if index == "1:mx":
            return "idx.x"
        if index == "1:nx":
            return "idx.x"

        if index == "1:my":
            return "idx.y"
        if index == "1:ny":
            return "idx.y"

        if index == "1:mz":
            return "idx.z"
        if index == "1:nz":
            return "idx.z"

        elif index == "l1:l2":
            return "idx.x"
        elif index == "l1+1:l2+1":
            return "idx.x+1"
        elif index == "l1+2:l2+2":
            return "idx.x+2"
        elif index == "l1+3:l2+3":
            return "idx.x+3"
        elif index == "l1-1:l2-1":
            return "idx.x-1"
        elif index == "l1-2:l2-2":
            return "idx.x-2"
        elif index == "l1-3:l2-3":
            return "idx.x-3"
        print("How to handle f index?: ",index)
        exit()
    def evaluate_integer(self,value):
        print("evaluating int->",value)
        #pencil specific values known beforehand
        value = value.replace("iux-iuu","0")
        value = value.replace("iuz-iuu","2")
        value = value.replace("l2-l1+1","nx")
        value = value.replace("n1","1+nghost")
        for mapping in self.flag_mappings:
            if mapping in value:
                value = replace_variable(value,mapping,self.flag_mappings[mapping])
        if value in self.known_values:
            orig_value = value
            while value in self.known_values:
                value = self.known_values[value]
            print("->",value)
            return value
        if "(" in value and ("/" in value or "*" in value):
            print("/ and * not supported with braces")
            exit()

        if "+" in value  and "*" not in value and "/" not in value:
            res = ""
            sum = 0
            for part in [part.strip() for part in value.split("+")]:
                if part.isnumeric() or (part[0] == "-" and part[:1].isnumeric()):
                    sum += eval(part)
                else:
                    res += part
            new_res = ""
            for part in [part.strip() for part in res.split("-")]:
                if part.isnumeric():
                    sum -= eval(part)
                else:
                    new_res += part
            new_res += "+" + str(sum)
            print("->",new_res)
            return new_res
        print("->",value)
        return value

    def evaluate_indexes(self,value):
        if value == ":":
            return ":"
        if ":" in value:
            lower,upper = [part.strip() for part in value.split(":")]
            lower = self.evaluate_integer(lower)
            upper = self.evaluate_integer(upper)
            return f"{lower}:{upper}"
        return self.evaluate_integer(value)

    def evaluate_boolean(self,value):
        if ".or." in value and ".and." in value:
            return value
        if value == "0./=0":
            return ".false."
        if value =="0.==0":
            return ".true."
        if value == "mx>nx":
            return ".true."
        if "(" in value:
            start_index = -1
            end_index = -1
            for i,char in enumerate(value):
                if char == "(" and start_index<0:
                    start_index = i
                if char == ")" and end_index<0:
                    end_index = i+1
            inside = self.evaluate_boolean(value[start_index+1:end_index-1])
            res_value = value[:start_index]
            res_value += inside
            res_value += value[end_index:]
            value = res_value
        if ".and." in value:
            parts = [self.evaluate_boolean(part.strip()) for part in value.split(".and.")]
            all_true = all([part == ".true." for part in parts])
            some_false = any([part == ".false." for part in parts])
            if all_true:
                return ".true."
            elif some_false:
                return  ".false."
            else:
                return value
        elif ".or." in value:
            parts = [self.evaluate_boolean(part.strip()) for part in value.split(".or.")]
            some_true= any([part == ".true." for part in parts])
            all_false= all([part == ".false." for part in parts])
            if some_true:
                return ".true."
            elif all_false:
                return  ".false."
            else:
                return value
        elif ".not." in value:
            res = self.evaluate_boolean(value.replace(".not.","").strip())
            if res == ".false.":
                return ".true."
            elif res == ".true.":
                return ".false."
            else:
                return res
        if "==" in value and len(value.split("==")) == 2:
            lhs, rhs = [part.strip() for part in value.split("==")]
            if lhs == rhs:
                return ".true."
            #integer comparison
            elif lhs.isnumeric() and rhs.isnumeric() and lhs != rhs:
                return ".false."
            elif lhs[0] == "'" and lhs[-1] == "'" and rhs[0] == "'" and rhs[-1] == "'":
                return ".false."
            #TOP /= BOT
            elif lhs == "top" and rhs == "bot":
                return ".false."
            elif lhs == "0." and rhs == "0.0":
                return ".true."
            else:
                return value
        if "/=" in value and len(value.split("/=")) == 2:
            lhs, rhs = [part.strip() for part in value.split("/=")]
            if lhs == rhs:
                return ".false."
            elif lhs.isnumeric() and rhs.isnumeric() and lhs != rhs:
                return ".true."
            elif lhs[0] == "'" and lhs[-1] == "'" and rhs[0] == "'" and rhs[-1] == "'":
                return ".true."
            elif lhs == "0." and rhs == "0.0":
                return ".false."
            else:
                #hack for checks if we are doing less than 3d runs which are not currently supported
                if lhs in ["nwgrid","nxgrid","nygrid","nzgrid"] and rhs in ["1"]:
                    return ".true."
                return value
        return value
    def find_module_files(self, module_name):
        return self.module_info[module_name]["files"]

    def find_subroutine_files(self, subroutine_name ):
        #Workaround since don't have time to copypaste since some files say they support functions which they don't support
        if "interface_funcs" in self.func_info[subroutine_name]:
            return self.func_info[subroutine_name]["files"]
        return [file for file in self.func_info[subroutine_name]["files"] if file in self.func_info[subroutine_name]["lines"]]

    def parse_module(self, module_name):
        if module_name not in self.parsed_modules:
            self.parsed_modules.append(module_name)
            file_paths = self.find_module_files(module_name)
            for file_path in file_paths:
                self.parse_file_for_static_variables(file_path)

    def get_chosen_modules(self,makefile):
        self.chosen_modules = []
        if makefile:
            lines = [line.strip().lower() for line in open(makefile,'r').readlines() if  line.strip() != "" and line.strip()[0] != "#" and line.split("=")[0].strip() != "REAL_PRECISION"] 
            for line in lines:
                if len(line.split("=")) == 2:
                    variable = line.split("=")[0].strip()
                    value = line.split("=")[1].strip()
                    self.chosen_modules.append(value)
                    if variable == "density":
                        self.chosen_modules.append(f"{value}_methods")

    def get_array_segments_in_line(self,line,variables):
        array_vars = self.get_arrays_in_line(line,variables)
        res = get_variable_segments(line, array_vars)
        return sorted(res, key = second_val)
    def get_struct_segments_in_line(self,line,variables):
        search_vars = []
        save_var = False
        buffer = ""
        for char in line:
            if char == "%":
                save_var = True
            if not(char.isalpha() or char.isnumeric()) and char != "%":
                if save_var:
                    search_vars.append(buffer.strip())
                    save_var = False
                buffer = ""
            else:
                buffer = buffer + char
        if save_var:
            search_vars.append(buffer.strip())
        return [seg for seg in get_variable_segments(line,search_vars) if seg[0] != ""]
    def get_arrays_in_line(self,line,variables):
        variables_in_line= get_used_variables_from_line(line)
        res = []
        for var in variables_in_line:
            if var in variables and len(variables[var]["dims"]) > 0:
                res.append(var)
        return res


    def parse_line(self, line):
        if len(line) == 0:
            return line
        if line[0] == "!":
            return line
        ## remove comment at end of the line
        iter_index = len(line)-1;
        num_of_single_quotes = 0
        num_of_double_quotes = 0
        for iter_index in range(len(line)):
            if line[iter_index] == "'":
                num_of_single_quotes += 1
            if line[iter_index] == '"':
                num_of_double_quotes += 1
            if line[iter_index] == "!" and num_of_single_quotes%2 == 0 and num_of_double_quotes%2 == 0 and (iter_index+1==len(line) or line[iter_index+1] != "="):
                line = line[0:iter_index]
                break                          
        return line.strip()
    def add_write(self, variable, line_num, is_local, filename, call_trace, line):
        self.static_writes.append({"variable": variable, "line_num": line_num, "local": is_local, "filename": filename, "call_trace": call_trace, "line": line})
    def get_pars_lines(self,prefix, name, lines):
        res = []
        in_pars = False
        for (line,count) in lines:
            line = line.strip().lower()
            if in_pars and line == "/":
                in_pars = False
            if in_pars:
                res.append(line)
            if f"&{name}{prefix}_pars" in line:
                in_pars = True
                res.append(line.replace(f"&{name}{prefix}_pars","").replace("/","").strip())
        #/ is replaced since it represents the end of line
        return [(line.replace("/","").strip(),0) for line in res]
    def get_lines(self, filepath, start=0, end=math.inf, include_comments=False):
        if filepath not in self.lines.keys():
            in_sub_name = None 
            lines = []
            start = False
            has_module_line = False
            has_program_line = False
            read_lines = open(filepath, 'r').readlines()
            index = 0
            while index<len(read_lines):
                start = False
                line = self.parse_line(read_lines[index].strip())
                match = re.search("include\s+(.+\.(h|inc))",line)
                if match and line[0] != "!":
                    header_filename = match.group(1).replace("'","").replace('"',"")
                    header_filepath = filepath.rsplit("/",1)[0].strip() + "/" + header_filename
                    if os.path.isfile(header_filepath):
                        header_lines = open(header_filepath, 'r').readlines()
                        read_lines = read_lines[:index+1] + header_lines + read_lines[1+index:]
                if len(line)>0:
                    #take care the case that the line continues after end of the line
                    if line[-1] == "&" and line[0] != "!":
                        while line[-1] == "&":
                            index += 1
                            next_line = read_lines[index].strip()
                            if len(next_line)>0 and next_line[0] != "!":
                                line = (line[:-1] + " " + self.parse_line(next_line)).strip()
                    # split multilines i.e. fsdfaf;fasdfsdf; into their own lines for easier parsing
                    for line in self.split_line(line):
                        if line[0] != '!':
                            search_line = line.lower().strip()

                            if "program" in search_line:
                                has_program_line = True
                            if search_line.split(" ")[0].strip() == "module" and search_line.split(" ")[1].strip() != "procedure":
                                has_module_line = True
                                module_name = search_line.split(" ")[1].strip()
                                if module_name not in self.module_info:
                                    self.module_info[module_name] = {}
                                if "files" not in self.module_info[module_name]:
                                    self.module_info[module_name]["files"] = []
                                self.module_info[module_name]["files"].append(filepath)
                                if filepath not in self.file_info:
                                    self.file_info[filepath] = {}
                                self.file_info[filepath]["module"] = module_name
                            if not in_sub_name:
                                search = re.search(f"\s?subroutine\s*([a-zA-Z0-9_-]*?)($|\s|\()", search_line)
                                if(search and "subroutine" in search_line):
                                    sub_name = search.groups()[0].strip()
                                    if sub_name not in self.func_info:
                                        self.func_info[sub_name] = {}
                                    if "lines" not in self.func_info[sub_name]:
                                        self.func_info[sub_name]["lines"] = {}
                                    self.func_info[sub_name]["lines"][filepath] = []
                                    if "files" not in self.func_info[sub_name]:
                                        self.func_info[sub_name]["files"] = []
                                    in_sub_name = sub_name
                                    start = True
                                search = re.search(f"\s?function\s*([a-zA-Z0-9_-]*?)($|\s|\()", search_line)
                                if(search and "function" in search_line):
                                    sub_name = search.groups()[0].strip()
                                    if sub_name not in self.func_info:
                                        self.func_info[sub_name] = {}
                                    if "lines" not in self.func_info[sub_name]:
                                        self.func_info[sub_name]["lines"] = {}
                                    self.func_info[sub_name]["lines"][filepath] = []
                                    if "files" not in self.func_info[sub_name]:
                                        self.func_info[sub_name]["files"] = []
                                    in_sub_name = sub_name
                                    start = True
                            if in_sub_name and not start:
                                if("end subroutine" in search_line or "endsubroutine" in search_line or "end function" in search_line or "endfunction" in search_line):
                                    self.func_info[in_sub_name]["lines"][filepath].append((self.parse_line(line),index))
                                    in_sub_name = None
                            if "interface" in line:
                                iter_index = index
                                res_line = line.lower()
                                if res_line.split(" ")[0].strip() == "interface" and len(res_line.split(" "))>1:
                                    sub_name = res_line.split(" ")[1].strip()
                                    if sub_name not in self.func_info:
                                        self.func_info[sub_name] = {}
                                    if "interface_funcs" not in self.func_info[sub_name]:
                                        self.func_info[sub_name]["interface_funcs"] = {}
                                    if "files" not in self.func_info[sub_name]:
                                        self.func_info[sub_name]["files"] = []
                                    self.func_info[sub_name]["files"].append(filepath)
                                    self.func_info[sub_name]["interface_funcs"][filepath] = []
                                    cur_index = index+1
                                    cur_line = read_lines[cur_index].lower()
                                    find = True
                                    while not ("end" in cur_line and "interface" in cur_line):
                                        if("module procedure" in cur_line):
                                            self.func_info[sub_name]["interface_funcs"][filepath].append(self.parse_line(cur_line.split("module procedure ")[1].strip()))
                                        elif("function" in cur_line):
                                            self.func_info[sub_name]["interface_funcs"][filepath].append(self.parse_line(cur_line.split("function")[1].split("(")[0].strip()))
                                        cur_index += 1
                                        cur_line = read_lines[cur_index].lower()
                        lines.append((line, index))
                        if in_sub_name and line[0] != '!':
                            self.func_info[in_sub_name]["lines"][filepath].append((self.parse_line(line),index))
                index += 1
            self.lines[filepath] = lines
            if filepath not in self.file_info:
                self.file_info[filepath] = {}
            self.file_info[filepath]["is_program_file"] = not has_module_line and has_program_line
        res = [x for x in self.lines[filepath] if x[1] >= start and x[1] <= end]
        if not include_comments:
            res = [x for x in res if x[0][0] != "!"] 
        return res


    def contains_subroutine(self, lines, function_name):
        for (line, count) in [(line[0].lower(),line[1]) for line in lines]:
            if "(" in function_name or ")" in function_name:
                return False
            if re.search(f"\s?subroutine {function_name}[\s(]",line) or re.search(f"\s?function {function_name}[\s(]",line) or re.search(f"interface {function_name}(\s|$|,)",line) or line==f"subroutine {function_name}" or line==f"function {function_name}":
                return True
        return False

    def get_used_modules(self,lines):
        modules = []
        for (line, count) in lines:
            if line.strip().split(" ")[0].strip() == "use":
                modules.append(line.strip().split(" ")[1].strip().split(",")[0].strip())
        return [module.lower() for module in modules]
    def get_own_module(self,filename):
        return self.file_info[filename]["module"]
        

    def get_mem_access_or_function_calls(self,lines):
        res = []
        for(line, count ) in filter(lambda x: not is_variable_line(x),lines):
            line = line.strip()
            matches = re.findall("[^'=' '\/+-.*()<>]+\(.+\)", line)
            if len(matches) > 0:
                res.extend(matches)
        return res


    def parse_variable(self, line_segment):
        iter_index = len(line_segment)-1
        end_index = iter_index
        start_index = 0
        num_of_left_brackets = 0
        num_of_right_brackets = 0
        while iter_index>0 and (line_segment[iter_index] in " ()" or num_of_left_brackets != num_of_right_brackets):
            elem = line_segment[iter_index]
            if elem == "(":
                num_of_left_brackets += 1
            elif elem == ")":
                num_of_right_brackets += 1
            iter_index -= 1
        end_index = iter_index+1

        while iter_index>0 and line_segment[iter_index] not in " *+-();":
            iter_index -= 1
        
        if iter_index == 0:
            res = line_segment[0:end_index]
        else:
            res = line_segment[iter_index+1:end_index]
        if "%" in res:
            end_index = 0
            while res[end_index] != "%":
                end_index += 1
            res = res[0:end_index]
        return res


    def get_writes_from_line(self,line,count,to_lower=True):
        res = []
        index = 0
        start_index = 0
        num_of_single_quotes = 0
        num_of_double_quotes = 0
        num_of_left_brackets = 0
        num_of_right_brackets = 0
        while index<len(line)-1:
            index += 1
            elem = line[index]
            if elem == "'":
                num_of_single_quotes += 1
            elif elem == '"':
                num_of_double_quotes += 1
            elif elem == "=" and line[index-1] not in "<>=!" and line[index+1] not in "<>=!" and num_of_single_quotes%2==0 and num_of_double_quotes%2==0 and num_of_left_brackets == num_of_right_brackets:
                write = self.parse_variable(line[start_index:index])
                variable = write
                if to_lower:
                    variable = variable.lower()
                res.append({"variable": variable.strip(), "line_num": count, "line": line, "is_static": write in self.static_variables, "value": line.split("=")[1].split(",")[0].strip()})
                start_index = index
            elif elem == "(":
                num_of_left_brackets += 1
            elif elem == ")":
                num_of_right_brackets += 1
        return res

    def get_writes(self,lines,exclude_variable_lines=True):
        res = []
        if exclude_variable_lines:
            lines = list(filter(lambda x: not is_variable_line(x), lines))
        for(line, count) in lines:
            res.extend(self.get_writes_from_line(line,count))
        return res
    def get_rhs_variable(self,line,to_lower=True):
        writes =  self.get_writes_from_line(line,0,to_lower)
        if len(writes) == 0:
            return None
        return writes[0]["variable"]
        # index = 0
        # start_index = 0
        # num_of_single_quotes = 0
        # num_of_double_quotes = 0
        # num_of_left_brackets = 0
        # num_of_right_brackets = 0
        # while index<len(line)-1:
        #     index += 1
        #     elem = line[index]
        #     if elem == "'":
        #         num_of_single_quotes += 1
        #     elif e:lem == '"':
        #         num_of_double_quotes += 1
        #     elif elem == "(":
        #         num_of_left_brackets = num_of_left_brackets + 1
        #     elif elem == ")":
        #         num_of_right_brackets = num_of_right_brackets + 1
        #     elif elem == "=" and line[index-1] not in "<>=!" and line[index+1] not in "<>=!" and num_of_single_quotes%2==0 and num_of_double_quotes%2==0 and num_of_left_brackets == num_of_right_brackets:

        #             return parse_variable(line[start_index:index])
        # return None
    def get_used_variables_from_line(self,line):
        characters_to_space= ["/",":", ",","+","-","*","(",")","=","<","!",">"]
        for character in characters_to_space:
            line = line.replace(character, " ")
        parts = [part.strip() for part in line.split(" ")]
        return [self.parse_variable(part) for part in parts if part]

    def get_used_variables(self,lines):
        res = []
        for(line, count) in filter(lambda x: not is_variable_line(x),lines):
            res.extend(self.get_used_variables_from_line(line))
        return res


    def get_function_name(self,line):
        iter_index = len(line)-1
        while(iter_index > 0):
            if line[iter_index] == "%":
                return line[0:iter_index]
            iter_index -= 1
        return line
    def get_function_calls_in_line(self, line, local_variables,exclude_variable_lines=True):
        return self.get_function_calls([(line,0)],local_variables, exclude_variable_lines)
    def get_function_calls(self,lines, local_variables,exclude_variable_lines=True):
        function_calls = []
        if exclude_variable_lines:
            lines = filter(lambda x: not is_variable_line(x),lines)
        for(line, count ) in lines:
            line = line.lower()
            function_call_indexes = []
            num_of_single_quotes = 0
            num_of_double_quotes = 0
            #get normal function calls i.e. with parameters and () brackets
            for i in range(len(line)):
                if line[i] == "'":
                    num_of_single_quotes += 1
                if line[i] == '"':
                    num_of_double_quotes += 1
                if line[i] == "(" and (i==0 or line[i-1] not in "^'=''\/+-.*\(\)<>^;:,") and num_of_double_quotes %2 == 0 and num_of_single_quotes %2 == 0:
                    function_call_indexes.append(i)
            for index in function_call_indexes:
                current_index = index-1
                if line[current_index] == " ":
                    while(line[current_index] == " "):
                        current_index -= 1
                while(line[current_index] not in " ^'=''\/+-.*\(\)<>^,;:" and current_index >= 0):
                    current_index -= 1
                current_index +=1
                function_name = self.get_function_name(line[current_index:index])
                save_index = current_index
                current_index = index
                function_name = function_name.strip()
                if function_name.lower() not in self.static_variables and function_name.lower() not in local_variables:
                    #step through (
                    current_index +=1
                    number_of_right_brackets = 1
                    number_of_left_brackets = 0
                    num_of_single_quotes = 0
                    num_of_double_quotes = 0
                    parameter_list_start_index = current_index
                    while(number_of_right_brackets>number_of_left_brackets):
                        parameter_list = line[parameter_list_start_index:current_index]
                        if line[current_index] == "'" and num_of_double_quotes %2 == 0:
                            num_of_single_quotes += 1
                        if line[current_index] == '"' and num_of_single_quotes %2 == 0:
                            num_of_double_quotes += 1
                        if line[current_index] == "(" and num_of_double_quotes %2 == 0 and num_of_single_quotes %2 == 0:
                            number_of_right_brackets += 1
                        elif line[current_index] == ")" and num_of_double_quotes %2 == 0 and num_of_single_quotes %2 == 0: 
                            number_of_left_brackets += 1
                        
                        current_index += 1
                     
                    parameter_list = line[parameter_list_start_index:current_index-1]
                    parameters = []
                    param ="" 
                    ##if inside brackets they are array indexes
                    num_of_left_brackets = 0
                    num_of_right_brackets= 0
                    num_of_single_quotes = 0
                    num_of_double_quotes = 0
                    for char in parameter_list:
                        if char == "'" and num_of_double_quotes %2 == 0:
                            num_of_single_quotes += 1
                        if char == '"' and num_of_single_quotes %2 == 0:
                            num_of_double_quotes += 1
                        if char == "(":
                            num_of_left_brackets = num_of_left_brackets + 1
                        if char == ")":
                            num_of_right_brackets = num_of_right_brackets + 1
                        if char == "," and num_of_left_brackets == num_of_right_brackets and num_of_double_quotes %2 == 0 and num_of_single_quotes %2 == 0: 
                            parameters.append(param.strip())
                            param = ""
                        else:
                            param = param + char
                    ## add last param
                    parameters.append(param.strip())
                    #if multithreading analysis
                    for i, param in enumerate(parameters):
                        if len(param) > 0 and param[0] == "-":
                            parameters[i] = param[1:]
                    function_name = function_name.strip()
                    if len(function_name) >0 and "%" not in function_name:
                        function_calls.append({"function_name": function_name, "parameters": [param.strip() for param in parameters if len(param.strip()) > 0], "range": (save_index,current_index), "line": line})
            
            #get function calls with call function i.e. call infront and no brackets
            function_call_indexes = []
            num_of_single_quotes = 0
            num_of_double_quotes = 0
            buffer = ""
            #get normal function calls i.e. with parameters and () brackets
            for i in range(len(line)):
                if line[i] == "'":
                    num_of_single_quotes += 1
                if line[i] == '"':
                    num_of_double_quotes += 1
                elif line[i] in "'()[];-+*/=^":
                    buffer = ""
                elif line[i] == " ":
                    if buffer == "call" and num_of_single_quotes %2 == 0 and num_of_double_quotes %2 == 0:
                        function_call_indexes.append(i)
                    buffer = ""
                else:
                    buffer = buffer + line[i]
            for index in function_call_indexes:
                #skip empty spaces
                while line[index] == " ":
                    index= index+1
                if line[index] == "=":
                    break
                start_index = index
                while index<len(line) and line[index] not in " (":
                    index = index+1
                function_name = line[start_index:index]
                if("=" in function_name):
                    print(line)
                    assert("=" not in function_name)
                #don't add duplicates, no need since no parameters
                if function_name not in [function_call["function_name"] for function_call in function_calls]:
                    function_name = function_name.strip()
                    ## % check is since the we can have for example p%inds(1) and %inds(1) would be marked as a function call.
                    if len(function_name) > 0 and "%" not in function_name:
                        function_calls.append({"function_name": function_name, "parameters": []})
            
         
        res =  [function_call for function_call in function_calls if not function_call["function_name"].isspace()]
        return res




    def get_contains_line_num(self, filename):
        lines = self.get_lines(filename)
        for (line, count) in lines:
            line = line.lower()
            if line.strip() == "contains":
                return count
            elif re.match("function\s+.+\(.+\)",line) or re.match("subroutine\s+.+\(.+\)",line):
                return count-1
        #whole file
        return lines[-1][1]

    def add_public_declarations_in_file(self,filename,lines):
        in_public_block = False
        for (line,count) in lines:
            if line == "public":
                in_public_block = True
            if line == "private":
                in_public_block = False
            parts = line.split("::")
            is_public_declaration = "public" == parts[0].split(",")[0].strip()
            if (is_public_declaration or in_public_block) and len(parts) == 2:
                variable_names = self.get_variable_names_from_line(parts[1])
                for variable in variable_names:
                    if variable in self.static_variables:
                        self.static_variables[variable]["public"] = True
            match = re.search("include\s+(.+\.h)",line)
            if match:
                header_filename = match.group(1).replace("'","").replace('"',"")
                directory = re.search('(.+)\/',filename).group(1)
                header_filepath = f"{directory}/{header_filename}"
                if os.path.isfile(header_filepath):
                    with open(header_filepath,"r") as file:
                        lines = file.readlines()
                        for line in lines:
                            parts = line.split("::")
                            is_public_declaration = "public" == parts[0].split(",")[0].strip()
                            if is_public_declaration:
                                variable_names = self.get_variable_names_from_line(parts[1])
                                for variable in variable_names:
                                    if variable in self.static_variables:
                                        self.static_variables[variable]["public"] = True

    def add_threadprivate_declarations_to_file(self,file, declaration):
        contents = []
        not_added = True
        for line in open(file, 'r').readlines():
            
            if (line.strip() == "contains" or "endmodule" in line) and not_added:
                not_added = False
                contents.append(f"{declaration}\n")
            elif not_added and re.match("function\s+.+\(.+\)",line) or re.match("subroutine\s+.+\(.+\)",line):
                contents.append(f"{declaration}\n")
                not_added = False
            contents.append(line)
        f = open(file, "w+")
        f.write("".join(contents))


    def get_variable_names_from_line(self,line_variables):
        res = []
        start_index=0
        end_index=0
        current_index=0
        parse_still = True
        while parse_still:
            current_elem = line_variables[current_index]
            if current_elem == "=":
            
                res.append(line_variables[start_index:current_index])
                parse_until_next_variable = True
                current_index +=1
                num_of_left_brackets = 0
                num_of_right_brackets = 0
                while parse_until_next_variable:
                    current_elem = line_variables[current_index]
                    not_inside_brackets = num_of_left_brackets == num_of_right_brackets
                    if current_elem == "!":
                        parse_until_next_variable = False
                        parse_still = False
                    if current_elem == "," and not_inside_brackets:
                        start_index = current_index+1
                        parse_until_next_variable=False
                    if current_elem == "(":
                        num_of_left_brackets += 1
                    if current_elem == ")":
                        num_of_right_brackets += 1
                    if parse_until_next_variable:
                        current_index +=1
                        if current_index >= len(line_variables):
                            start_index = current_index+1
                            parse_until_next_variable = False
                            parse_still = False
            elif current_elem == ",":
                res.append(line_variables[start_index:current_index])
                start_index = current_index + 1
            elif current_elem == "!":
                parse_still = False
                res.append(line_variables[start_index:current_index+1])
            current_index += 1
            if current_index >= len(line_variables):
                parse_still= False
                if current_index > start_index:
                    res.append(line_variables[start_index:current_index+1])
        return [x.replace("&","").replace("!","").replace("(:)","").strip() for x in res if x.replace("&","").replace("!","").strip() != ""]
    def parse_file_for_struct(self,file_path,struct):
        lines = self.get_lines(file_path)
        in_type = False 
        for (line,count) in lines:
            line = line.lower().strip()
            if line == "contains":
                in_type=False
                break
            if f"type {struct}" in line and "(" not in line and ")" not in line and "end" not in line:
                in_type = True
                if struct not in self.struct_table:
                    self.struct_table[struct] = {}
            elif in_type and ("endtype" in line or "end type" in line):
                in_type = False
                break
            elif in_type:
                self.get_variables([(line,count)], self.struct_table[struct], file_path)
    def parse_struct(self, struct):
        for module_name in self.parsed_modules:
            file_paths = self.find_module_files(module_name)
            for file_path in file_paths:
                self.parse_file_for_struct(file_path,struct)

    def get_variables(self, lines, variables, filename):
        module = None
        if "module" in self.file_info[filename]:
            module = self.get_own_module(filename)
            if module not in self.module_variables:
                self.module_variables[module] = {}
        in_type = False
        for i, (line, count) in enumerate(lines):
            parts = [part.strip() for part in line.split("::")]
            is_variable = True
            start =  parts[0].strip()
            type = start.split(",")[0].strip()
            if type == "public":
                is_variable = False
                if len(parts) == 2:
                    variable_names = [name.lower() for name in self.get_variable_names_from_line(parts[1])]
                    for variable in variable_names:
                        if variable in self.func_info:
                            if filename not in self.func_info[variable]["files"]:
                                self.func_info[variable]["files"].append(filename)
                                #in case is interface func add the interface funcs it has in the file
                                if "lines" not in self.func_info[variable]:
                                    for interface_func in self.func_info[variable]["interface_funcs"][filename]:
                                        if filename not in self.func_info[interface_func]["files"]:
                                            self.func_info[interface_func].append(filename)
            if "end" not in line and (line.split(" ")[0].split(",")[0].strip() == "type" and len(line.split("::")) == 1) or (line.split(" ")[0].split(",")[0].strip() == "type" and "(" not in line.split("::")[0]):
                in_type = True
            if in_type and ("endtype" in line.lower().strip()  or "end type" in line.lower().strip()):
                in_type = False
            if is_variable and not in_type and len(parts)>1:
                allocatable = "allocatable" in [x.strip() for x in start.split(",")]
                saved_variable = "save" in [x.strip() for x in start.split(",")]
                public = "public" in [x.strip() for x in start.split(",")]
                is_parameter = "parameter" in [x.strip() for x in start.split(",")]
                is_optional = "optional" in [x.strip() for x in start.split(",")]
                dimension = []
                if "dimension" in start:
                    _, end = [(m.start(0), m.end(0)) for m in re.finditer("dimension", start)][0]
                    while start[end] != "(":
                        end = end + 1
                    rest_of_line = start[end+1:]
                    num_of_left_brackets = 1
                    num_of_right_brackets= 0
                    index = ""
                    i = 0
                    while num_of_left_brackets != num_of_right_brackets:
                        char = rest_of_line[i]
                        if char == "(":
                            num_of_left_brackets = num_of_left_brackets + 1
                        if char == ")":
                            num_of_right_brackets = num_of_right_brackets + 1
                        if char == "," and num_of_left_brackets == num_of_right_brackets + 1: 
                            dimension.append(index.strip())
                            index = ""
                        else:
                            index = index + char
                        i = i +1
                    index = index[:-1]
                    dimension.append(index.strip())
                line_variables = parts[1].strip()
                variable_names = [name.lower() for name in self.get_variable_names_from_line(line_variables)]
                for i, variable_name in enumerate(variable_names):
                    dims = dimension
                    if "(" in variable_name:
                        variable_name = variable_name.split("(")[0].strip()
                        search = re.search(f"{variable_name}\(((.*?))\)",line) 
                        ## check if line is only specifying intent(in) or intent(out)
                        if search :
                            dims = [index.strip() for index in search.group(1).split(",")]
                    ## get first part of .split(" ") if f.e. character len(something)
                    type = type.split(" ")[0].strip()
                    #filter intent(in) and intent(inout)
                    if "intent(" not in type and ".true." not in variable_name and ".false." not in variable_name:
                        if variable_name not in variables:
                            #if struct type parse more to get the type
                            if "type" in type:
                                struct_type = re.search("\((.+?)\)",line).group(1)
                                type = struct_type
                                if struct_type not in self.struct_table:
                                    self.parse_struct(struct_type.lower())
                            if "(" in type:
                                type = type.split("(")[0].strip()
                            if "rkind8" in line and type == "real":
                                type = "double"
                            type = type.lower()
                            variables[variable_name] = {"type": type, "dims": dims, "allocatable": allocatable, "origin": [filename], "public": public, "threadprivate": False, "saved_variable": (saved_variable or "=" in line_variables.split(",")[i]), "parameter": is_parameter, "on_target": False, "optional": is_optional}
                        else:
                            variables[variable_name]["origin"].append(filename)
                            variables[variable_name]["origin"] = list(set(variables[variable_name]["origin"]))
                        if module is not None:
                            if variable_name not in self.module_variables[module]:
                                if "type" in type:
                                    struct_type = re.search("\((.+?)\)",line).group(1)
                                    type = struct_type
                                    if struct_type not in self.struct_table:
                                        self.parse_struct(struct_type.lower())
                                if "(" in type:
                                    type = type.split("(")[0].strip()
                                if "rkind8" in line and type == "real":
                                    type = "double"
                                type = type.lower()
                                if variable_name == "t" and type != "double":
                                    #TODO fix this
                                    print("HMM","t line")
                                    print(line)
                                    print(filename)
                                    print(in_type)
                                else:
                                    self.module_variables[module][variable_name] = {"type": type, "dims": dims, "allocatable": allocatable, "origin": [filename], "public": public, "threadprivate": False, "saved_variable": (saved_variable or "=" in line_variables.split(",")[i]), "parameter": is_parameter, "on_target": False, "optional": is_optional}
                            else:
                                self.module_variables[module][variable_name]["origin"].append(filename)
                                self.module_variables[module][variable_name]["origin"] = list(set(variables[variable_name]["origin"]))


                            

        return variables




    def parse_file_for_static_variables(self, filepath):
        modules = self.get_always_used_modules(filepath)
        for module in modules:
            self.parse_module(module)
        self.load_static_variables(filepath)

    def add_threadprivate_declarations_in_file(self,filename):
        res = []
        lines = self.get_lines(filename, include_comments=True)
        index = 0
        while index<len(lines):
            line, _ = lines[index]
            if "!$omp" in line and "threadprivate" in line.lower():
                if line[-1] == "&":
                    if line[-1] == "&":
                        while line[-1] == "&":
                            index += 1
                            next_line = lines[index][0].strip().replace("!$omp","")
                            if len(next_line)>0:
                                line = (line[:-1] + " " + next_line).strip()
                search = re.search("(threadprivate|THREADPRIVATE)\((.*)\)",line)
                if search:
                    variable_names = [variable.strip() for variable in search.group(2).split(",")]
                    for variable in variable_names:
                        res.append(variable)
                        if variable in self.static_variables:
                            self.static_variables[variable]["threadprivate"] = True
            index = index+1
        return res
    
    def function_is_already_declared(self,function,filename):
        res = []
        lines = self.get_lines(filename, include_comments=True)
        index = 0
        while index<len(lines):
            line, _ = lines[index]
            if "!$omp" in line and "declare target" in line.lower():
                if line[-1] == "&":
                    if line[-1] == "&":
                        while line[-1] == "&":
                            index += 1
                            next_line = lines[index][0].strip().replace("!$omp","")
                            if len(next_line)>0:
                                line = (line[:-1] + " " + next_line).strip()
                search = re.search("(declare target)\((.*)\)",line)
                if search:
                    function_names = [variable.strip() for variable in search.group(2).split(",")]
                    if function in function_names:
                        return True
            index = index+1
        return False
    def add_declare_target_declarations_in_file(self,filename):
        res = []
        lines = self.get_lines(filename, include_comments=True)
        index = 0
        while index<len(lines):
            line, _ = lines[index]
            if "!$omp" in line and "declare target" in line.lower():
                if line[-1] == "&":
                    if line[-1] == "&":
                        while line[-1] == "&":
                            index += 1
                            next_line = lines[index][0].strip().replace("!$omp","")
                            if len(next_line)>0:
                                line = (line[:-1] + " " + next_line).strip()
                search = re.search("(declare target)\((.*)\)",line)
                if search:
                    variable_names = [variable.strip() for variable in search.group(2).split(",")]
                    for variable in variable_names:
                        res.append(variable)
                        for var in self.static_variables:
                            if var == variable or f"{var}_offload" == variable:
                                self.static_variables[var]["on_target"] = True

            index = index+1
        return res


    def load_static_variables(self, filename, dst=None):
        if dst is None:
            dst = self.static_variables
        self.parsed_files_for_static_variables.append(filename)
        static_variables_end = self.get_contains_line_num(filename)
        self.get_variables(self.get_lines(filename, 1, static_variables_end), dst, filename)
        self.add_public_declarations_in_file(filename,self.get_lines(filename, 1, static_variables_end))
        self.add_threadprivate_declarations_in_file(filename)
        self.add_declare_target_declarations_in_file(filename)
        return dst

    def get_subroutine_variables(self, filename, subroutine_name):
        return self.get_variables(self.get_subroutine_lines(subroutine_name,filename))

    def get_always_used_modules(self, filename):
        module_declaration_end = self.get_contains_line_num(filename)
        return [module_name for module_name in self.get_used_modules(self.get_lines(filename, 1, module_declaration_end)) if module_name.lower() not in self.ignored_modules]

    def get_subroutine_modules(self, filename, subroutine_name):
        return [module_name for module_name in  self.get_used_modules(self.get_subroutine_lines(subroutine_name,filename)) if module_name in self.ignored_modules]

    def get_all_modules_in_file_scope(self,filename,modules):
        modules_to_add = self.get_always_used_modules(filename)
        modules_to_add.append(self.get_own_module(filename))
        already_added_modules = modules.copy()
        modules.extend(modules_to_add)
        for module in modules_to_add:
            if module not in already_added_modules:
                for file in self.find_module_files(module):
                    self.get_all_modules_in_file_scope(file, modules)
        return list(set(modules))
    def get_all_modules_in_subroutine_scope(self,filename,subroutine_name):
        if filename not in self.modules_in_scope:
            self.modules_in_scope["filename"] = self.get_all_modules_in_file_scope(filename, [self.get_own_module(filename)])
        res = list(set(self.modules_in_scope["filename"] + self.get_subroutine_modules(filename,subroutine_name)))
        res = [module for module in res if module not in self.ignored_modules]
        for module in res:
            self.parse_module(module)
        return res

    
    def get_local_module_variables(self,filename,subroutine_name):
        res = {}
        #Variables in own file take precedence 
        self.load_static_variables(filename,res)

        for module in self.get_all_modules_in_subroutine_scope(filename,subroutine_name):
            if module not in self.module_variables and module in self.parsed_modules:
                print("module was parsed but not self.module_variables",module)
                exit()
            for var in self.module_variables[module]:
                if var not in res:
                    res[var] = self.module_variables[module][var]
        return res

    def get_subroutine_lines(self, filename, subroutine_name):
        func_info = self.get_function_info(subroutine_name)
        if "lines" not in func_info:
            print("trying to get func lines without lines", subroutine_name, filename)
            exit()
        return func_info["lines"][filename]

    def get_parameters(self,line):
        line = line.lower()
        check_subroutine= re.search("\s?subroutine.+\((.+)\)",line)
        if check_subroutine:
            res =  [parameter.split("=")[-1].split("(")[0].strip().lower() for parameter in check_subroutine.group(1).split(",")]
            return res
        else:
            ##The return parameter needs to be removed
            if "result" in line:
                line = line.split("result")[0].strip()

            param_search = re.search(".?function.+\((.+)\)",line)
            if not param_search:
                return []
            res =  [parameter.split("=")[-1].split("(")[0].strip().lower() for parameter in param_search.group(1).split(",")]
            return [res_part.lower() for res_part in res]

    def get_parameter_mapping(self, parameters, parameter_list):
        # return [x[0] for x in enumerate(parameter_list)]
        # Commented out for testing
        mapping = []
        for i, param in enumerate(parameter_list):
            if param[4] is None:
                mapping.append(i)
            else:
                for j,sub_param in enumerate(parameters):
                    if sub_param == param[4]:
                        mapping.append(j)
        return mapping

    def get_static_parameters(self,line,parameter_list):
        parameters = self.get_parameters(line)
        res = []
        mapping = self.get_parameter_mapping(parameters,parameter_list)
        already_done = []
        for i, map_index in enumerate(mapping):
            res.append((parameters[map_index], parameter_list[i][:-1]))
            already_done.append(map_index)
        for i in range(len(parameters)):
            if i not in already_done:
                res.append((parameters[i],("",False,"",[])))
        return res
    def get_var_info_from_array_access(self,parameter,local_variables,local_module_variables):
        var = parameter[0].split("(")[0].strip()
        if "%" in var:
            var_name, field = var.split("%",1)
            if var_name in local_variables:
                sizes = self.struct_table[local_variables[var_name]["type"]][field]["dims"]
            else:
                sizes = self.struct_table[self.static_variables[var_name]["type"]][field]["dims"]
        else:
            if var in local_variables:
                sizes = local_variables[var]["dims"]
            else:
                sizes = self.static_variables[var]["dims"]
        if parameter[0] == "(":
            parameter[0] == parameter[0][1:]
        if parameter[-1] == ")":
            parameter[0] == parameter[0][:-1]
        indexes = get_indexes(parameter[0],var,0)
        ##count the number of looped indexes
        dim = 0
        dims = []
        for i, index in enumerate(indexes):
            ##can have inline array as array range
            if ":" in index:
                if index == ":":
                    dims.append(sizes[i])
                else:
                    lower,upper = [part.strip() for part in index.split(":")]
                    dims.append(f"{upper}-{lower}+1")
            elif index.replace("(","").replace(")","").strip()[0] == '/' and index.replace("(","").replace(")","").strip()[-1] == '/':
                dims.append(":")
        source =  local_variables if var in local_variables else local_module_variables
        is_static = source[var]["saved_variable"] or var in local_module_variables if "%" not in var else False
        # type =  res[2] if res[2] != "" else source[var]["type"]
        type = source[var]["type"] if "%" not in var else "pencil_case"
        return (var,is_static, type, dims)
    def get_param_info(self,parameter,local_variables,local_module_variables):
        if len(parameter[0]) == 0:
            print("INCORRECT PARAM",parameter)
            exit()
        if parameter[0] in local_variables:
            return (parameter[0],parameter[1],local_variables[parameter[0]]["type"],local_variables[parameter[0]]["dims"])
        if parameter[0] in local_module_variables:
            return (parameter[0],parameter[1],local_module_variables[parameter[0]]["type"],local_module_variables[parameter[0]]["dims"])
        if parameter[0] in self.static_variables:
            return (parameter[0],parameter[1],self.static_variables[parameter[0]]["type"],self.static_variables[parameter[0]]["dims"])
        if parameter[0][0] == "(":
            parameter = (parameter[0][1:],parameter[1])
        if parameter[0] == ".true." or parameter[0] == ".false.":
            return (parameter[0],False,"logical",[])
        is_sum = False
        is_product = False
        is_division = False
        is_difference = False
        num_of_left_brackets = 0
        num_of_right_brackets= 0
        possible_var = ""
        in_array = False
        for char_index,char in enumerate(parameter[0]):
            if char == "(":
                num_of_left_brackets += 1
                if parameter[0][char_index-1] == " ":
                    back_index = char_index-1
                    while(parameter[0][back_index] == " " and back_index>0):
                        back_index -= 1
                    if back_index != 0 and parameter[0][back_index] not in "*-+/(":
                        in_array = True
                else:
                    if char_index != 0 and parameter[0][char_index-1] not in "*-+-/(":
                        # in_array = possible_var in local_variables or possible_var in local_module_variables 
                        in_array = True
                possible_var = ""
            elif char == ")":
                num_of_right_brackets += 1
                if num_of_left_brackets == num_of_right_brackets:
                    in_array = False
                possible_var = ""
            elif char == "+" and not in_array:
                is_sum = True
                possible_var = ""
            elif char == "*" and not in_array:
                is_product= True
                possible_var = ""
            elif char == "-" and not in_array:
                is_difference= True
                possible_var = ""
            elif char == "/" and not in_array and parameter[0][char_index-1] != "(" and parameter[0][char_index+1] != ")" and char_index != 0 and char_index != len(parameter[0])-1:
                if char_index < len(parameter[0])-1:
                    is_division = parameter[0][char_index+1] not in ")!"
                else:
                    is_division= True
                possible_var = ""
            else:
                possible_var = possible_var + char
        operations = (is_sum,is_difference,is_product,is_division)
        print("PARAM",parameter,operations)
        #Inline array
        if not is_division and parameter[0].replace("(","")[0] == "/" and parameter[0].replace(")","")[-1] == "/":
            par_str = parameter[0].strip()
            if par_str[0] != "(":
                par_str = "(" + par_str
            if par_str[1] == "/":
                par_str = par_str[0] + par_str[2:]
            if par_str[-1] != ")":
                par_str = par_str + ")"
            if par_str[-2] == "/":
                par_str = par_str[:-2] + par_str[-1]
            par_str = "inline_array" + par_str
            parameters = self.get_function_calls_in_line(par_str, local_variables)[0]["parameters"]
            info = self.get_param_info((parameters[0],False),local_variables,local_module_variables)
            return (parameter[0],"False",info[2],[":"])
        func_calls = self.get_function_calls_in_line(parameter[0],local_variables)
        if len(func_calls)>0 and not any(operations):
            first_call = func_calls[0]
                #Functions that simply keep the type of their arguments
            if first_call["function_name"] in ["sqrt","alog","log","exp","sin","cos","log","abs"]:
                return self.get_param_info((first_call["parameters"][0],False),local_variables,local_module_variables)
            #Array Functions that return single value if single param else an array
            if first_call["function_name"] in ["sum"]:
                new_param = (first_call["parameters"][0],False)
                inside_info =  self.get_param_info(new_param,local_variables,local_module_variables)
                if len(first_call["parameters"]) == 1:
                    return (first_call["parameters"][0],False,inside_info[2],[])
                else:
                    return (first_call["parameters"][0],False,inside_info[2],[":"])
            #Array Functions that return scalar value, multiple params, return type is passed on first param
            if first_call["function_name"] in ["dot_product"]:
                
                inside_info =  self.get_param_info((first_call["parameters"][0],False),local_variables, local_module_variables)
                return (parameter[0],False,inside_info[2],[])
            #Array Functions that return the largest value in params, multiple params, return type is passed on first param
            if first_call["function_name"] in ["max","min"]:
                inside_info =  self.get_param_info((first_call["parameters"][0],False),local_variables, local_module_variables)
                return (parameter[0],False,inside_info[2],inside_info[3])
            
            if first_call["function_name"] in ["maxval","minval"]:
                inside_info =  self.get_param_info((first_call["parameters"][0],False),local_variables, local_module_variables)
                return (parameter[0],False,inside_info[2],inside_info[3][:-1])

            #SPREAD
            if first_call["function_name"] == "spread":
                first_param = parameter[0].split("(",1)[1].split(")")[0].split(",")[0].strip()
                inside_info =  self.get_param_info((first_param,False),local_variables, local_module_variables)
                spread_dims = [":" for x in inside_info[3]]
                #Result dim is source array dim+1
                spread_dims.append(":")
                return (parameter[0],False,inside_info[2],spread_dims)
            if first_call["function_name"] == "real":
                inside_info = self.get_param_info((first_call["parameters"][0],False),local_variables, local_module_variables)
                return (parameter[0],False,"real",inside_info[3])
            if first_call["function_name"] == "int":
                inside_info = self.get_param_info((first_call["parameters"][0],False),local_variables, local_module_variables)
                return (parameter[0],False,"integer",inside_info[3])
            if first_call["function_name"] == "trim":
                return (parameter[0],False,"character",[])
            if first_call["function_name"] == "len":
                return (parameter[0],False,"integer",[])
            print("unsupported func")
            print(first_call)
            exit()
        elif parameter[0] in local_variables:
            return (parameter[0],parameter[1],local_variables[parameter[0]]["type"],local_variables[parameter[0]]["dims"])
        elif parameter[0] in local_module_variables:
            return (parameter[0],parameter[1],local_module_variables[parameter[0]]["type"],local_module_variables[parameter[0]]["dims"])
        elif "'" in parameter[0] or '"' in parameter[0]:
            return (parameter[0],parameter[1],"character",[])
        elif any(operations):
            op_chars = []
            if is_sum:
                 op_chars.append("+")
            if is_difference:
                 op_chars.append("-")
            if is_product:
                 op_chars.append("*")
            if is_division:
                 op_chars.append("/")
            parts = [part.strip() for part in split_by_ops(parameter[0].strip(),op_chars)]
            # print("PARTS",parts)
            parts_res = ("",False,"",[])
            for part in parts:
                if len(part) > 0  and part[0] == "(" and sum([char == "(" for char in part]) > sum([char == ")" for char in part]):
                    part = part[1:]
                if len(part) > 0  and part[-1] == ")" and sum([char == "(" for char in part]) < sum([char == ")" for char in part]):
                    part = part[:-1]
                part_res = self.get_param_info((part,False),local_variables,local_module_variables)
                if parts_res[2] == "" or len(part_res[3])>len(parts_res[3]):
                    parts_res = (parts_res[0],False,part_res[2],part_res[3])
            return parts_res
        elif "%" in parameter[0] and "'" not in parameter[0] and '"' not in parameter[0]:
            var_name,field_name = [part.strip() for part in parameter[0].split("%")]
            ##var_name can be array access if array of structures
            var_name = var_name.split("(")[0]
            struct = ""
            if var_name in local_variables:
                struct = local_variables[var_name]["type"]
            else:
                struct = local_module_variables[var_name]["type"]
            field_name = field_name
            print("field name",field_name)
            if "(" in field_name: 
                var_dims = self.get_var_info_from_array_access(parameter,local_variables,local_module_variables)[-1]
                print("FIELD DIMS",var_dims)
            else:
                field = self.struct_table[struct][field_name]
                var_dims = field["dims"]
            field_name = field_name.split("(")[0]
            field = self.struct_table[struct][field_name]
            return (var_name,parameter[1],field["type"],var_dims)
        elif ".and." in parameter[0] or ".not." in parameter[0] or ".or." in parameter[0]:
            return (parameter[0], False, "logical", [])
        elif "(" in parameter[0] and ")" in parameter[0] and not any(operations) and ".and." not in parameter[0] and ".not." not in parameter[0]:
            var = parameter[0].split("(")[0].strip()
            if var in local_variables or var in local_module_variables:
                return self.get_var_info_from_array_access(parameter,local_variables,local_module_variables)
            ##Boolean intrinsic funcs
            elif var in ["present","isnan","associated","allocated","all","any"]:
                return (parameter[0],False,"logical",[])
            elif var in ["trim"]:
                return (parameter[0],False,"character",[])
            else:
                #check if function in source code
                file_paths = self.find_subroutine_files(var)
                if len(file_paths)>0:
                    interfaced_functions = self.get_interfaced_functions(file_paths[0],var)
                    if len(interfaced_functions)>0:
                        sub_lines = self.get_subroutine_lines(var,file_paths[0])
                        declaration_line = sub_lines[0]
                        if "logical" in declaration_line[0]:
                            return (parameter[0],False,"logical",[])
                        ##inline array
                        elif "real" in declaration_line[0]:
                            res_type = "real"
                        elif "(/" in parameter[0]:
                            type = "real" if "." in parameter[0].split("/")[1].split(",")[0] else "integer"
                            return (parameter[0],False,type,[":"])
                        else:
                            print(f"DIDN'T FIND RETURN TYPE OF {var} :(")
                            print(parameter[0])
                            exit()
                else:
                    print(f"DIDN'T FIND VAR {var} :(")
                    print(parameter[0])
                    print(is_sum, is_product,is_division,is_difference)
                    exit()
        else:
            type =""
            if "'" in parameter[0] or '"' in parameter[0]:
                type = "character"

            elif parameter[0].isnumeric() or parameter[0][1:].isnumeric():
                type = "integer"
            elif parameter[0].replace(".","").isnumeric() or parameter[0][1:].replace(".","").isnumeric():
                type = "real"
            elif "." in parameter:
                if all([part.isnumeric() or part[1:].isnumeric() for part in parameter[0].split(",")]):
                    type = "real"
                #is scientific number
                if parameter[0].replace(".","").replace("-","").replace("e","").replace("+","").isnumeric():
                    return (parameter[0],False,"real",[])
            return (parameter[0],parameter[1],type,[])
    def get_static_passed_parameters(self,parameters,local_variables,local_module_variables):
        original_parameters = parameters
        parameters = [parameter.lower() for parameter in parameters]
        for i,param in enumerate(parameters):
            if len(param)>0 and param[0] == "-":
                parameters[i] = param[1:]
        res = list(zip(parameters,list(map(lambda x: (x in local_module_variables and x not in local_variables) or (x in local_variables and local_variables[x]["saved_variable"]),parameters))))
        for i,parameter in enumerate(res):
            param_name = parameter[0]
            func_calls = self.get_function_calls_in_line(param_name,local_variables)
            if len(func_calls) == 0:
               param_name = parameter[0].split("=")[-1]
            info = self.get_param_info((param_name,parameter[1]),local_variables,local_module_variables)
            if len(parameter[0].split("=")) == 2:
                param_name,passed_param = [part.strip().lower() for part in parameter[0].split("=")]
                res[i] = (info[0],info[1],info[2],info[3],param_name)
            else:
                res[i] = (info[0],info[1],info[2],info[3],None)

        return res

    def get_interfaced_functions(self,file_path,subroutine_name):
        func_info = self.get_function_info(subroutine_name)
        if("interface_funcs" not in func_info):
            return [subroutine_name.lower()]
        if(file_path not in func_info["interface_funcs"]):
            #colliding names for a subroutine and an interface
            return [subroutine_name.lower()]
        if(len(func_info["interface_funcs"][file_path]) == 0):
             return [subroutine_name.lower()]
        res = [sub.lower() for sub in func_info["interface_funcs"][file_path]]
        if subroutine_name in res:
            print("subroutine in it's own interface")
            print(subroutine_name)
            exit()
        return [sub for sub in res if file_path in self.func_info[sub]["lines"]]

    def generate_save_array_store(self,store_variable):
        res = f"{store_variable}_generated_array(imn"
        if len(self.static_variables[store_variable]['dims']) == 0:
            res += ",1"
        else:
            for i in range(len(self.static_variables[store_variable]['dims'])):
                res += ",:"
        res += f") = {store_variable}\n"
        return res

    def generate_read_from_save_array(self,store_variable):
        res = f"{store_variable} = {store_variable}_generated_array(imn"
        if len(self.static_variables[store_variable]['dims']) == 0:
            res += ",1"
        else:
            for i in range(len(self.static_variables[store_variable]['dims'])):
                res += ",:"
        res +=")\n"
        return res

    def generate_allocation_for_save_array(self,store_variable):
        res = f"{self.static_variables[store_variable]['type']}, dimension ("
        if self.static_variables[store_variable]["allocatable"]:
            res += ":"
        else:
            res += "nx*ny"
        if len(self.static_variables[store_variable]['dims']) == 0:
            if self.static_variables[store_variable]["allocatable"]:
                res += ",:"
            else:
                res += ",1"
        else:
            for dimension in self.static_variables[store_variable]["dims"]:
                res += f",{dimension}"
        res += ")"
        if self.static_variables[store_variable]["allocatable"]:
            res += ", allocatable"
        res += f" :: {store_variable}_generated_array\n"
        return res

    def save_static_variables(self):
        with open("static_variables.csv","w",newline='') as csvfile:
            writer = csv.writer(csvfile)
            for variable in self.static_variables.keys():
                writer.writerow([variable,self.static_variables[variable]["type"],self.static_variables[variable]["dims"],self.static_variables[variable]["allocatable"],self.static_variables[variable]["origin"],self.static_variables[variable]["public"],self.static_variables[variable]["threadprivate"]])
        
    def read_static_variables(self):
        with open("static_variables.csv","r",newline='') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                self.static_variables[row[0]] = {"type": row[1], "dims": [dim for dim in row[2].replace("'","").strip('][').split(', ') if dim != ""], "allocatable": (row[3].lower() in ("yes", "true", "t", "1")), "origin": row[4], "public": (row[5].lower() in ("yes", "true", "t", "1")), "threadprivate": (row[6].lower() in ("yes", "true", "t", "1"))}

    def save_static_writes(self,static_writes):
        keys = static_writes[0].keys()
        with open("writes.csv","w",newline='') as output_file:
            dict_writer = csv.DictWriter(output_file, keys)
            dict_writer.writeheader()
            dict_writer.writerows(static_writes)

    def read_static_writes(self):
        with open("writes.csv", mode="r") as infile:
            return [dict for dict in csv.DictReader(infile)]

    def get_threadprivate_declarations_and_generate_threadpriv_modules(self,files,threadprivate_variables,critical_variables):
        modules_file= open(f"{self.directory}/threadpriv_modules.inc","w")
        use_modules_file= open(f"{self.directory}/use_threadpriv_modules.inc","w")
        copy_in_call_file = open(f"{self.directory}/copyin.inc","w")
        threadprivate_declarations = {}
        exceptions = ["cdata.f90","cparam.f90"]
        for file in files:
            threadprivate_variables_to_add = []
            static_variables_end = self.get_contains_line_num(file)
            vars = list(self.get_variables(self.get_lines(file, 1, static_variables_end), {}, file).keys())
            # print("vars")
            # for var in vars:
            #    print(var)
                # print(self.static_variables[var])
                # print(self.static_variables[var]["threadprivate"])
            vars_to_make_private = [variable for variable in vars if variable in self.static_variables and (variable in threadprivate_variables or  self.static_variables[variable]["threadprivate"]) and variable not in critical_variables and variable != ""]
            if(len(vars_to_make_private) > 0):
                module = self.get_own_module(file)
                copy_in_file = open(f"{file.replace('.f90','')}_copyin.inc","w")
                if not any([exp in file for exp in exceptions]):
                    copy_in_file.write(f"subroutine copyin_{module}\n")
                copy_in_file.write("!$omp parallel copyin( &\n")
                var_num = 0
                for var in vars_to_make_private:
                    if var_num > 0:
                        copy_in_file.write(f"!$omp ,{var}  &\n")
                    else:
                        copy_in_file.write(f"!$omp {var}  &\n")
                    var_num += 1
                res_lines = [f"!$omp threadprivate({var})" for var in vars_to_make_private]
                res_line = ""
                for x in res_lines:
                    res_line += x + "\n"
                threadprivate_declarations[file] = res_line
                modules_file.write(f"{self.get_own_module(file)}\n")
                use_modules_file.write(f"use {self.get_own_module(file)}\n")
        
                copy_in_file.write("!$omp )\n")
                copy_in_file.write("!$omp end parallel\n")
                if not any([exp in file for exp in exceptions]):
                    copy_in_file.write(f"end subroutine copyin_{module}\n")
                    copy_in_call_file.write(f"call copyin_{module}\n")
                copy_in_file.close()

        modules_file.close()
        copy_in_call_file.close()
        use_modules_file.close()
        return threadprivate_declarations

    def generate_copy_in(self,files,threadprivate_variables):
        variables = []
        for file in files:
            threadprivate_variables_to_add = []
            static_variables_end = self.get_contains_line_num(file)
            vars = list(self.get_variables(self.get_lines(file, 1, static_variables_end), {}, file).keys())
            private_vars = [variable for variable in vars if variable in self.static_variables and (variable in threadprivate_variables or self.static_variables[variable]["threadprivate"]) and variable != ""]
            variables.extend(private_vars)
            variables = list(set(variables))
        print("Creating copyin.inc")
        split_str = ",\n!$omp"
        copyin_line = f"!$omp parallel copyin({split_str.join(variables)})\n"
        copyin_file = open("copyin.inc", "w")
        # copyin_file.write("subroutine copyin_func()\n")
        copyin_file.write(copyin_line)
        copyin_file.write("!$omp end parallel")
        copyin_file.close()

                    
    def add_atomic_declarations(self,variables):
        files = self.used_files
        for file in files:
            res_contents = []
            contents = self.get_lines(file,include_comments=True)
            for i,(line,count) in enumerate(contents):
                res_contents.append(line)
                if line[0] != "!":
                    writes = self.get_writes_from_line(line,count)

                    #handle where writes
                    if "if" in line or line.split("(")[0].strip() == "where":
                        for variable in variables:
                            if variable in [write["variable"] for write in writes]:
                                last_line = res_contents[-1]
                                res_contents[-1] = "!omp critical\n"
                                res_contents.append(last_line)
                                res_contents.append("!omp end critical\n")
                    #handle do loop writes
                    elif re.match("do\s+.=.+,\s?",line.split(";")[0]) and len(line.split(";"))==2:
                        for variable in variables:
                            if variable in [write["variable"] for write in writes]:
                                res_contents[-1] = (f"{line.split(';')[0]}\n")
                                res_contents.append("!omp critical\n")
                                res_contents.append(f"{line.split(';')[1]}\n")
                                res_contents.append("!omp end critical\n")
                                #let's see if there is a corresponding enddo
                                current_index = i+1
                                no_end_do = True
                                iter_line = contents[current_index][0]
                                while no_end_do and not re.match("do\s+.=.+,\s?",iter_line):
                                    no_end_do = not (re.match("enddo",iter_line) or re.match("end do",iter_line))
                                    current_index += 1
                                    iter_line = contents[current_index][0]
                                if no_end_do:
                                    res_contents.append("enddo")
                                
                    else:
                        for variable in variables:
                            if variable in [write["variable"] for write in writes]:
                                last_line = res_contents[-1]
                                res_contents[-1] = "!omp critical\n"
                                res_contents.append(last_line)
                                res_contents.append("!omp end critical\n")
                                #If one one's the more performant omp atomic
                                #last_line = res_contents[-1]
                                #res_contents[-1] = "!$omp atomic\n"
                                #res_contents.append(last_line)

            with open(f"./out/{file}","w") as f:
                f.write("\n".join(res_contents))

    def make_variables_public(self,variables):
        files = list(set([self.static_variables[variable]["origin"] for variable in variables ]))
        for file in files:
            res_contents = []
            contents = self.get_lines(file,include_comments=True)
            in_module_declaration = True
            in_type = False
            for (line,count) in contents:
                res_contents.append(line)
                if line.split(" ")[0] == "type":
                    in_type = True
                if line.split(" ")[0] == "endtype":
                    in_type = False
                if in_module_declaration and not in_type:
                    if re.match("function\s+.+\(.+\)",line) or re.match("subroutine\s+.+\(.+\)",line) or line.strip() == "contains":
                        in_module_declaration = False
                    parts = line.split("::")
                    if len(parts) > 1:
                        line_variables = parts[1].strip()
                        variable_names = self.get_variable_names_from_line(line_variables)
                        res_contents.extend([f"public :: {variable}" for variable in variable_names if variable in variables and self.static_variables[variable]["origin"] == file and not self.static_variables[variable]["public"]])

            with open(f"./out/{file}","w") as f:
                f.write("\n".join(res_contents))
    def get_function_info(self,function):
        return self.func_info[function]
    def is_elemental(function):
        func_info = self.get_function_info(function)
        return func_info["is_elemental"]

    def choose_correct_interfaced_function(self,subroutine_name,interfaced_functions,parameter_list,file_path):
        #Unfortunately different modules have different calling structures of the same function so can be 0 for a file.
        suitable_functions = []
        for function in interfaced_functions:
            is_suitable = True
            subroutine_lines = self.get_subroutine_lines(function, file_path)
            parameters = self.get_parameters(subroutine_lines[0][0])
            is_elemental = "elemental" in subroutine_lines[0][0]
            local_variables = {parameter:v for parameter,v in self.get_variables(subroutine_lines, {},file_path).items() }
            mapping = self.get_parameter_mapping(parameters,parameter_list)
            print("FUNCTION",function, "IN", file_path)
            print("PARAMS")
            for param in parameters:
                print(param, local_variables[param])
            #if mapping is less than parameter list than some named optional paramaters are not present in sub parameters
            is_suitable  = len(mapping) == len(parameter_list) and len(parameters) >= len(parameter_list)
            ##check if type and length of dims match between passed parameter and function parameter.
            ## All other parameters need to be optional 
            if is_suitable:
                print(len(parameter_list))
                for i, passed_param in enumerate(parameter_list):
                    is_suitable =  is_suitable and  passed_param[2] == local_variables[parameters[mapping[i]]]["type"] 
                    if not is_elemental:
                        is_suitable = is_suitable and len(passed_param[3]) == len(local_variables[parameters[mapping[i]]]["dims"])
                for i in [j for j in range(len(parameters)) if j not in mapping]:
                    is_suitable = is_suitable and local_variables[parameters[i]]["optional"]
            if is_suitable:
                num_of_needed_params = 0
                for param in parameters:
                    if not local_variables[param]["optional"]:
                        num_of_needed_params = num_of_needed_params + 1
                is_suitable = is_suitable and len(parameter_list) >= num_of_needed_params
            if is_suitable:
                suitable_functions.append(function)
        num_of_suitable_needed = 1 if self.offloading else 0
        if len(suitable_functions) > 1 or len(suitable_functions)<num_of_suitable_needed:
            print(f"There are {len(suitable_functions)} suitable functions for the interface call: ",subroutine_name, "in file", file_path)
            print("Params: ",parameter_list)
            print("Original candidates: ", interfaced_functions)
            print("Possibilities", suitable_functions)
            exit()
        return suitable_functions
    def choose_right_module(self,filepaths):
        print("CHOOSING RIGHT MODULE")
        for i, path in enumerate(filepaths):
            for module in self.chosen_modules:
                if path.lower() == f"{self.directory}/{module}.f90":
                    return filepaths[i]
        ##For testing return only if in Makefile file
        ##no module is default if not found module
        # for path in filepaths:
        #     if f"{self.directory}/no" in path:
        #         return path
        print("did not found module in files",filepaths)
        exit()
        
    def parse_subroutine_all_files(self, subroutine_name, call_trace, check_functions, offload,local_variables,file_called_from, layer_depth=math.inf, parameter_list=[], only_static=True):
        if layer_depth<0:
            return []
        self.subroutine_order += 1
        if subroutine_name not in self.subroutine_modifies_param:
            self.subroutine_modifies_param[subroutine_name] = {}
        file_paths = self.find_subroutine_files(subroutine_name)
        print("Parse all files", subroutine_name, parameter_list)
        #if no files than it must be in the file_called_from
        if len(file_paths) == 0:
            file_paths = [file_called_from]
        global_modified_list = []
        param_types = [(param[2],len(param[3])) for param in parameter_list]
        all_functions = []
        all_found_functions = []
        for file_path in file_paths:
            interfaced_functions = self.get_interfaced_functions(file_path,subroutine_name)
            all_functions.extend(interfaced_functions)
            interfaced_functions = self.choose_correct_interfaced_function(subroutine_name,interfaced_functions,parameter_list,file_path)
            for function in interfaced_functions:
                all_found_functions.append(function)
                self.parse_subroutine_in_file(file_path, function, check_functions,offload,global_modified_list, layer_depth, call_trace,parameter_list,only_static)
        if len(all_found_functions) == 0 and len(all_functions)>0:
            print(f"There are no suitable function for the interface call: ",subroutine_name)
            print("Params: ",parameter_list)
            print("Original candidates: ", all_functions)
            print(self.static_variables["t"])
            exit()
        self.subroutine_modifies_param[subroutine_name][str(param_types)] = global_modified_list
    def replace_var_in_lines(self, lines, old_var, new_var):
        # return [add_splits(replace_variable(line,old_var,new_var)) for line in lines]
        return [replace_variable(line,old_var,new_var) for line in lines]
    def get_subroutine_lines(self,subroutine_name,filename):
        func_info = self.get_function_info(subroutine_name)
        print("sub name",subroutine_name)
        print("files", [file for file in func_info["lines"]])
        print("files", [file for file in func_info["files"]])
        print("lines", [file for file in func_info["files"]])
        return func_info["lines"][filename]
    def turn_f_size(self,params):
        if len(params) > 1:
            dim = params[1][0]
            if dim == "1":
                return "mx"
            elif dim == "2":
                return "my"
            elif dim  == "3":
                return "mz"
            else:
                print("Weird dim in size")
                print(line)
                exit()
        else:
            return "mx*my*mz"
    def get_dim_info(self,param,local_variables,variables_in_scope,writes):
        if param in local_variables:
            src = local_variables
        else:
            src = self.static_variables
        if not any([":" in dim or "size" in dim for dim in src[param]["dims"]]):
            print("here",param)
            return src[param]["dims"]
        var_writes = []
        res_dims = [":" for i in range(len(src[param]["dims"]))]
        for write in writes:
            if write["variable"] == param:
                var_writes.append(write)
                vars_in_line = [var for var in get_used_variables_from_line(write["line"]) if var in self.static_variables or var in local_variables]
                vars_dims = []
                for var in vars_in_line:
                    if var in local_variables:
                        src = local_variables
                    else:
                        src = self.static_variables
                    vars_dims.append(src[var]["dims"])
                dims_index = 0
                for var_dims in vars_dims:
                    for i, var_dim in enumerate(var_dims):
                        if var_dim != ":" and "size" not in var_dim:
                            if len(res_dims) <= i:
                                res_dims.append(var_dim)
                            else:
                                res_dims[i] = var_dim
        print("returning",res_dims)
        return res_dims

    def get_size(self, func_call, local_variables, variables_in_scope, writes):
        params = func_call["parameters"]
        if params[0] == "f":
            return self.turn_f_size(params)
        elif params[0] in local_variables and all([":" not in dim and "size" not in dim for dim in local_variables[params[0]]['dims']]):
            sizes = local_variables[params[0]]["dims"]
        elif params[0] in local_variables and len(local_variables[params[0]]["dims"]) == 1 and "size" in local_variables[params[0]]["dims"][0]:
            size_func_call = self.get_function_calls_in_line(local_variables[params[0]]["dims"][0],local_variables)[0]
            sizes = self.get_size(size_func_call, local_variables, variables_in_scope, writes)
            local_variables[params[0]]["dims"] == ["*".join(sizes)]
        else:
            print("size call", func_call)
            info = self.get_param_info((params[0],False),local_variables,self.static_variables)
            print("param info",info)
            sizes = info[3]
        if len(func_call["parameters"]) == 1:
            return "*".join(sizes)
        else:
            dim = func_call["parameters"][1]
            return  sizes[int(dim)-1]
    def replace_func_call(self, line,func_call, replacement):
        return line[:func_call["range"][0]] + replacement + line[func_call["range"][1]:]
    def get_replaced_body(self, original_subroutine_name, filename, parameter_list, function_call_to_replace, variables_in_scope,global_init_lines):
        ##in case is interfaced call get the correct subroutine
        print("GETTING REPLACED BODY FOR: ", function_call_to_replace)
        interfaced_functions = self.get_interfaced_functions(filename,original_subroutine_name)
        if len(interfaced_functions)>1:
            interfaced_functions = self.choose_correct_interfaced_function(original_subroutine_name,interfaced_functions,parameter_list,filename)
        subroutine_name = interfaced_functions[0]
        if "inlined_lines" not in self.func_info[subroutine_name]:
            self.inline_all_function_calls(filename,subroutine_name)
        if filename not in self.func_info[subroutine_name]["inlined_lines"]:
            self.inline_all_function_calls(filename,subroutine_name)
        subroutine_lines = self.func_info[subroutine_name]["inlined_lines"][filename]
        #subroutine_lines = self.get_subroutine_lines(subroutine_name,filename)

        print("FILENAME: ",filename,original_subroutine_name,subroutine_name)
        init_lines = [line[0] for line in subroutine_lines if is_init_line(line[0])]

        for module in self.get_used_modules(subroutine_lines):
            self.parse_module(module)
            self.parsed_modules.append(module)
            self.parsed_modules = list(set(self.parsed_modules))
        lines = subroutine_lines
        is_function = "function" in subroutine_lines[0][0]
        type = None
        local_variables = {parameter:v for parameter,v in self.get_variables(subroutine_lines, {},filename).items() }
        result_func_calls = [call for call in self.get_function_calls_in_line(subroutine_lines[0][0],local_variables) if call["function_name"] == "result"]
        has_named_return_value= is_function and "result" in subroutine_lines[0][0] and len(result_func_calls) == 1
        if has_named_return_value:
            return_var = result_func_calls[0]["parameters"][0]
        elif is_function:
            #if not named them default name is the function name
            return_var = subroutine_name

        if is_function:
            if "real" in subroutine_lines[0][0]:
                type = "real"
            elif "logical" in subroutine_lines[0][0]:
                type = "logical"
        #return value type is the local variable with it's name
        if is_function and not type:
            if return_var in local_variables:
                type = local_variables[return_var]["type"]
            else:
                if "character" in subroutine_lines[0][0]:
                    type = "character"
                elif "real" in subroutine_lines[0][0]:
                    type = "real"
                elif "integer" in subroutine_lines[0][0]:
                    type = "integer"
        # function_calls = self.get_function_calls(subroutine_lines[1:], local_variables)
        # for func_call in function_calls:
        #     if func_call["function_name"] not in self.ignored_subroutines:
        #         func_call_index = get_function_call_index(func_call["function_name"],lines)
        #         replaced_body_lines = self.expand_function_call(lines, subroutine_name, filename, func_call["function_name"])
        #         local_res_lines = lines[:func_call_index] + replaced_body_lines + lines[func_call_index+1:]
        #         lines = local_res_lines
        res_init_lines = []
        for line in init_lines:
            if len([part.strip() for part in line.split("::")[-1].split(",")]) > 1:
                res = line.split("::")[-1].strip()
                vars = []
                num_of_left_brackets = 0
                num_of_right_brackets = 0
                buffer  = ""
                for char in res:
                    if char == "(":
                        num_of_left_brackets += 1
                    if char == ")":
                        num_of_right_brackets += 1
                    if char == "," and num_of_right_brackets == num_of_left_brackets:
                        vars.append(buffer.strip())
                        buffer = ""
                    else:
                        buffer = buffer + char
                if buffer.strip() != "":
                    vars.append(buffer.strip())
                for var in vars:
                    res_init_lines.append(line.split("::")[0].strip() + "::" + var)
            else:
                res_init_lines.append(line)
        init_lines = res_init_lines
        params = self.get_parameters(subroutine_lines[0][0])
        remove_indexes = []
        for i,line in enumerate(init_lines):
            for param in params:
                if param in [part.strip() for part in line.split("::")[-1].split(",")]:
                    remove_indexes.append(i)
        init_lines = [x[1] for x in enumerate(init_lines) if x[0] not in remove_indexes]
        print("init lines after",init_lines) 
        print([[part.strip() for part in line.split("::")[-1].split(",")] for line in init_lines])

        new_lines = [line[0].lower() for line in lines]
        mapping = self.get_parameter_mapping(params, parameter_list)

        passed_param_names = [x.split("=")[-1].strip() for x in function_call_to_replace["parameters"]]
        present_params = [x[1] for x in enumerate(params) if x[0] in mapping]
        not_present_params = [param for param in params if param not in present_params]

        #remove_lines that write to parameters that are not passed
        remove_indexes = []
        for i,line in enumerate(new_lines):
            writes = self.get_writes_from_line(line,0)
            if(len(writes) == 1):
                if writes[0]["variable"] in not_present_params:
                    remove_indexes.append(i)

        new_lines = [x[1] for x in enumerate(new_lines) if x[0] not in remove_indexes]
        #replace variables with passed values
        for i, passed_param in enumerate(function_call_to_replace["parameters"]):
            #in case was passed as a named param
            passed_param = passed_param.split("=")[-1].strip()
            init_lines = self.replace_var_in_lines(init_lines, params[mapping[i]], passed_param)
            new_lines = self.replace_var_in_lines(new_lines, params[mapping[i]], passed_param)

        #replace all local_params that are not passed params with function specific names to prevent name collision with other inlined funcs
        for local_variable in [x for x in local_variables if x not in present_params]:
            init_lines = self.replace_var_in_lines(init_lines, local_variable, f"{local_variable}_{self.inline_num}")
            new_lines = self.replace_var_in_lines(new_lines, local_variable,f"{local_variable}_{self.inline_num}")

        if is_function:
            new_lines = self.replace_var_in_lines(new_lines,return_var,f"{return_var}_return_value_{self.inline_num}")
        if "step" in subroutine_name:
            print("new lines")
            print(subroutine_name, " is_function?: ", is_function)
            print(subroutine_name, " type?: ", type)
            for line in new_lines:
                print(line)
        #replace present with false if not present and true if present
        for i,line in enumerate(new_lines):
            present_func_calls = [func_call for func_call in self.get_function_calls_in_line(line,local_variables) if func_call["function_name"] == "present"]
            while len(present_func_calls) > 0:
                func_call = present_func_calls[0]
                present_var = func_call["parameters"][0]
                line = self.replace_func_call(line, func_call, ".true." if present_var in [par[0] for par in parameter_list] else ".false.")
                present_func_calls = [func_call for func_call in self.get_function_calls_in_line(line,local_variables) if func_call["function_name"] == "present"]
            new_lines[i] = line
        #replace loptest with false if not present
        for i,line in enumerate(new_lines):
            loptest_func_calls = [func_call for func_call in self.get_function_calls_in_line(line,local_variables) if func_call["function_name"] == "loptest"]
            while len(loptest_func_calls) > 0:
                print("in loptest loop")
                func_call = loptest_func_calls[0]
                present_var = func_call["parameters"][0]
                if len(func_call["parameters"]) == 1:
                    if present_var not in [par[0] for par in parameter_list]:
                        print("replacing loptest")
                        print("before ",line)
                        line = self.replace_func_call(line, func_call, ".false.")
                        print("after",line)
                    else:
                        print("replacing loptest")
                        print("before ",line)
                        line = self.replace_func_call(line, func_call, present_var)
                        print("after",line)
                loptest_func_calls= [func_call for func_call in self.get_function_calls_in_line(line,local_variables) if func_call["function_name"] == "loptest"]
            new_lines[i] = line

        init_variables= {parameter:v for parameter,v in self.get_variables([(line,0) for line in init_lines], {},filename).items() }
        for init_var in init_variables:
            if init_var == "tmp":
                print("what went wrong got tmp without underscore")
                print(function_call_to_replace)
                print(init_lines)
                exit()
        global_init_lines.extend(init_lines)
        global_init_lines = list(set(global_init_lines))
        lines = [(line,0) for line in new_lines]
        #Todo has to evaluate whether return line is hit or not
        has_return_line = False 
        for (line,count) in lines:
            if line.strip() == "return":
                has_return_line = True
        if has_return_line:
            print("\n\nlines before elimination\n\n")
            for (line,count) in lines:
                print(line)
            lines = [(line,0) for line in self.eliminate_while([x[0] for x in lines])]
        if has_return_line:
            print("\n\nlines after elimination\n\n")
            for (line,count) in lines:
                print(line)
            if self.has_no_ifs(lines):
                remove_indexes = []
                has_return_line = False
                for line_index,(line,count) in enumerate(lines):
                    has_return_line = has_return_line or line == "return"
                    if has_return_line:
                        remove_indexes.append(line_index)
                lines = [x[1] for x in enumerate(lines) if x[0] not in remove_indexes]
            else:
                print("don't know what to do since return statement in a conditional branch")
                exit()
        for (line,count) in lines:
            if "if(idiag_u2mr==0) return" in line:
                print("don't want to inline this func")
                print(subroutine_name, original_subroutine_name)
                exit()
        return ([line for line in lines if is_body_line(line[0])],is_function,type)
    def parse_subroutine_in_file(self, filename, subroutine_name, check_functions, offload, global_modified_list = [], layer_depth=math.inf, call_trace="", parameter_list=[], only_static=True):
        print("parse_in_file", subroutine_name, filename,parameter_list)
        if layer_depth < 0:
            return []
        if subroutine_name not in self.subroutine_modifies_param:
            self.subroutine_modifies_param[subroutine_name] = {}
        subroutine_lines = self.get_subroutine_lines(subroutine_name, filename)
        lines = subroutine_lines[1:]
        own_module = self.get_own_module(filename)
        self.parse_module(self.get_own_module(filename))
        for module in self.get_subroutine_modules(filename, subroutine_name):
            self.parse_module(module)
        local_variables = {parameter:v for parameter,v in self.get_variables(lines, {},filename).items() }
        local_module_variables = self.get_local_module_variables(filename,subroutine_name)
        print("param_list",parameter_list)
        parameters = self.get_static_parameters(subroutine_lines[0][0], parameter_list)
        if self.subroutine_order == 0:
            call_trace = subroutine_name
        else:
            call_trace = f"{call_trace} -> {subroutine_name}"
        writes = self.get_writes(lines)
        param_is_modified= {}
        for (parameter,(passed_parameter,is_static,_,_)) in parameters:
            param_is_modified[parameter] = {"is modified": False, "line" : "", "filename": filename, "call trace": call_trace}
        for write in writes:
            ##Remember to check for function return values
            if (write["variable"] not in self.static_variables) and (write["variable"].strip() not in local_variables) and write["variable"].strip() != subroutine_name.strip():
                print("Variable: ", write["variable"], "In function: ", subroutine_name, "Didn't have an origin")
                print([var for var in  local_variables])
                print([var for var in  self.static_variables])
                print("Either a bug or a missing file!")
                exit()
            write["is_static"] = (write["variable"] in self.static_variables and write["variable"] not in local_variables) or (write["variable"] in local_variables and local_variables[write["variable"]]["saved_variable"])
            for (parameter,(passed_parameter,is_static,_,_)) in parameters:
                if write["variable"] == parameter:
                    write["variable"] = passed_parameter
                    write["is_static"] = is_static
                    param_is_modified[parameter] = {"is modified": True, "line" : write["line"], "filename": filename, "call trace": call_trace}
        for i, (param,(passed_parameter,is_static,_,_)) in enumerate(parameters):
            if i < len(global_modified_list):
                if param_is_modified[param]["is modified"]:
                    global_modified_list[i] = param_is_modified[param]
            else:
                global_modified_list.append(param_is_modified[param])
        for write in writes: 
            if write["is_static"]:
                self.add_write(write["variable"], write["line_num"], write["variable"] in local_variables, filename, call_trace, write["line"])
        for function_call in self.get_function_calls(lines, local_variables):
            function_name = function_call["function_name"].lower()
            parse = True
            if function_name in check_functions or function_name.lower().startswith("mpi"):
                self.found_function_calls.append((call_trace, function_name, parameter_list))
                parse = False
            parse = parse and function_name.lower().strip() not in self.ignored_subroutines
            if parse:
                print("FINDING PARAMS IN", subroutine_name, filename, "FOR", function_name)
                print(f"FUNC PARAMS for {function_name}",function_call["parameters"])
                new_param_list = self.get_static_passed_parameters(function_call["parameters"],local_variables,local_module_variables)
                param_types = [(param[2],len(param[3])) for param in new_param_list]
                parse = (parse and ((function_name, str(param_types)) not in self.parsed_subroutines) and function_name != subroutine_name)
        # for i in range(len(new_param_list)):
            if parse:
                self.parsed_subroutines.append((function_name,str(param_types)))
                self.parse_subroutine_all_files(function_name,call_trace, check_functions, offload, local_variables, filename, layer_depth-1,new_param_list,only_static)
                #see if passed params were modified
                #if modified add them as writes and in case they were themselves passed params signal that they are modified for subroutines calling this subroutine
                for i in range(len(new_param_list)):
                    for j, (parameter,(passed_parameter,is_static,_,_)) in enumerate(parameters):
                        if new_param_list[i][0] == parameter:
                            if self.subroutine_modifies_param[function_name][str(param_types)][i]["is modified"]:
                                global_modified_list[j] = {"is modified": True, "filename": self.subroutine_modifies_param[function_name][str(param_types)][i]["filename"], "line": self.subroutine_modifies_param[function_name][str(param_types)][i]["line"], "call trace": self.subroutine_modifies_param[function_name][str(param_types)][i]["call trace"]}
                                if is_static:
                                    self.add_write(passed_parameter, 0, False, filename, self.subroutine_modifies_param[function_name][str(param_types)][i]["call trace"],self.subroutine_modifies_param[function_name][str(param_types)][i]["line"])
            elif function_name.lower().strip() not in self.ignored_subroutines and function_name not in check_functions and not function_name.lower().startswith("mpi"):

                # check if the subroutine name is itself in case of an recursive call
                if function_name == subroutine_name:
                    for i in range(len(new_param_list)):
                        if global_modified_list[i]["is modified"] and new_param_list[i][0] in self.static_variables or (new_param_list[i][0] in local_variables and local_variables[new_param_list[i][0]]["saved_variable"]):
                            self.add_write(new_param_list[i][0], 0, new_param_list[i][0] in local_variables, filename, call_trace, global_modified_list[i]["line"])
                        for j, (parameter,(passed_param,is_static,_,_)) in enumerate(parameters):
                            if new_param_list[i][0] == parameter and global_modified_list[i]["is modified"]:
                                global_modified_list[j] = {"is modified": True, "filename": global_modified_list[i]["filename"], "line": global_modified_list[i]["line"], "call trace": call_trace}
                                if is_static:
                                    self.add_write(passed_parameter, 0, False, filename, call_trace, global_modified_list[i]["line"])
                else:
                    new_param_list = self.get_static_passed_parameters(function_call["parameters"],local_variables,local_module_variables)
                    for i in range(len(new_param_list)):
                        if self.subroutine_modifies_param[function_name][str(param_types)][i]["is modified"] and new_param_list[i][0] in self.static_variables or (new_param_list[i][0] in local_variables and local_variables[new_param_list[i][0]]["saved_variable"]): 
                            self.add_write(new_param_list[i][0], 0, new_param_list[i][0] in local_variables, self.subroutine_modifies_param[function_name][str(param_types)][i]["filename"], self.subroutine_modifies_param[function_name][str(param_types)][i]["call trace"], self.subroutine_modifies_param[function_name][str(param_types)][i]["line"])
    def evaluate_ifs(self,lines,local_variables):
        for line_index, line in enumerate(lines):
            # print("evaluating ifs",line)
            orig_line = line
            check_line = line.replace("then","").replace("else if","if",1).replace("elseif","if",1).replace("if","call if",1)
            has_then = "then" in line
            if_func_calls = list(filter(lambda func_call: func_call["function_name"] == "if", self.get_function_calls_in_line(check_line,local_variables)))
            replacement_value = None
            if len(if_func_calls) == 1:
                func_call = if_func_calls[0]
                if len(func_call["parameters"]) == 1 and func_call["parameters"][0] in local_variables and "value" in local_variables[func_call["parameters"][0]]:
                    if not replacement_value and local_variables[func_call["parameters"][0]]["value"] == ".false.":
                        replacement_value = ".false."
                    if not replacement_value and local_variables[func_call["parameters"][0]]["value"] == ".true.":
                        replacement_value = ".true."
                if len(func_call["parameters"]) == 1 and not replacement_value:
                    for variable in [variable for variable in local_variables if "value" in local_variables[variable]]:
                        if variable in func_call["parameters"][0]:
                            if not replacement_value:
                                replacement_value = func_call["parameters"][0]
                            replacement_value = replace_variable(replacement_value, variable, local_variables[variable]["value"])
                if not replacement_value:
                    replacement_value = self.evaluate_boolean(func_call["parameters"][0])
                replacement_value = self.evaluate_boolean(replacement_value)
                #     if not replacement_value and ".and." in func_call["parameters"][0]:
                        # parts = [part.strip() for part in func_call["parameters"][0].split()]
                        # if all([part == ".true." for part in parts]):
                        #     replacement_value = ".true."
                        # if any([part == ".false." for part in parts]):
                        #     replacement_value = ".false."
                if replacement_value:
                    
                    replacement_value = self.evaluate_boolean(replacement_value)
                    if "else" not in line:
                        if_str = "if"
                    else:
                        if_str ="else if"
                    orig_line = line 
                    lines[line_index] = self.replace_func_call(check_line, func_call, f"{if_str}({replacement_value})").replace("call if","if").replace("call else if", "else if")
                    if "if(f /) then" in lines[line_index]:
                        print("I fucked up in evaluate index")
                        exit()
                    if "elseelse if" in lines[line_index]:
                        print("Fdsfsdafsd")
                        print(orig_line)
                        print(check_line)
                        print(if_str)
                        exit()
                    if has_then:
                        lines[line_index] = lines[line_index] + "then"
                        if "if(f /) then" in lines[line_index]:
                            print(lines[line_index])
                            print(orig_line)
                            print("I fucked up in evaluate index")
                            exit()
    def global_loop_replacer(self,segment,segment_index,line,local_variables,info):
        if segment[0] in local_variables:
            dims = local_variables[segment[0]]["dims"]
        else:
            dims = self.static_variables[segment[0]]["dims"]
        indexes = get_segment_indexes(segment,line,len(dims))
        if info["iterator"] in indexes:
            for i, index in enumerate(indexes):
                if info["iterator"] in index:
                    # if i!=info["replacement_index"]:
                    #     print("wrong replacement_index")
                    #     print(line)
                    #     print(info)
                    #     exit()
                    new_lower= index.replace(info["iterator"],info["replacement_lower"])
                    new_upper= index.replace(info["iterator"],info["replacement_upper"])
                    if new_lower == "1" and new_upper == dims[i]:
                        indexes[i] = ":"
                    else:
                        indexes[i] = new_lower + ":" + new_upper
        return build_new_access(segment[0],indexes)
    def remove_global_loops(self, lines,local_variables,global_loop_lines,iterators):
        remove_indexes = []
        variables = merge_dictionaries(local_variables,self.static_variables)
        for loop_arr in global_loop_lines:
            remove_indexes.append(loop_arr[0][1])
            remove_indexes.append(loop_arr[-1][1])

            #this needed because funny business with nsave
            iterator_index = loop_arr[0][1] - 1
            rhs_var = self.get_rhs_variable(lines[iterator_index])
            while(rhs_var == "n" or rhs_var == "m"):
                remove_indexes.append(iterator_index)
                iterator_index -= 1
                rhs_var = self.get_rhs_variable(lines[iterator_index])
            iterator_index = loop_arr[-1][1] + 1
            rhs_var = self.get_rhs_variable(lines[iterator_index])
            while(rhs_var == "n" or rhs_var == "m"):
                remove_indexes.append(iterator_index)
                iterator_index += 1
                rhs_var = self.get_rhs_variable(lines[iterator_index])


        for loop_arr_index, loop_arr in enumerate(global_loop_lines):
            print(iterators)
            iterator,replacement_lower, replacement_upper, replacement_index = iterators[loop_arr_index]
            for line,line_index in loop_arr:
                lines[line_index] = self.replace_segments(self.get_array_segments_in_line(line,variables),line,self.global_loop_replacer,local_variables,{
                    "iterator": iterator, 
                    "replacement_lower": replacement_lower,
                    "replacement_upper": replacement_upper,
                    "replacement_index": replacement_index,
                })
        return [x[1] for x in enumerate(lines) if x[0] not in remove_indexes]
    def is_vector(self, lines, vector, local_variables):
        if vector in local_variables:
            src = local_variables
        elif vector in self.static_variables:
            src = self.static_variables
        else:
            print("I have to ask")
            print(vector)
            exit()
            return False 
        if src[vector]["dims"] not in [["nx","3"]]:
            return False
        return True
    def inline_replacer(self,segment,segment_index,line,local_variables,info):
        if segment_index== 0:
            return line[segment[1]:segment[2]]
        else:
            var = segment[0].strip()
            if var in local_variables and len(local_variables[var]["dims"]) <= 1 and var in info["possible_values"]:
                indexes = get_segment_indexes(segment,line,1)
                is_safe = indexes[0] in [":","l1:l2"]
                if is_safe:
                    res = "("  + replace_variable(line[segment[1]:segment[2]].split("(",1)[0].strip(), var, info["possible_values"][var]) + ")"
                    return res
                    # add_line =  replace_variable(line[segment[1]:segment[2]], var, possible_values[var])
                else:
                    print("not safe abort")
                    print(line)
                    print(line[segment[1]:segment[2]])
                    exit()
            res = line[segment[1]:segment[2]]
            return res
    def replace_segments(self,segments,line,map_func,local_variables,info):
        res_line = ""
        last_index = 0 
        for sg_index, segment in enumerate(segments):
            res_line = res_line + line[last_index:segment[1]]
            res_line = res_line + map_func(segment,sg_index,line,local_variables,info) 
            last_index = segment[2]
        res_line = res_line + line[last_index:]
        return res_line
    def inline_1d_writes(self, lines,local_variables):
            analyse_lines = self.get_analyse_lines(lines)
            writes = self.get_writes([(line[0],line[3]) for line in analyse_lines])
            # vars_to_check_safety = ["tmp1","tmp2"]
            vars_to_check_safety = local_variables
            safe_vars_to_inline = []
            for var in vars_to_check_safety:
                var_writes= [write for write in writes if write["variable"] == var]

                no_unsafe_writes = False 
                # for write in var_writes:
                #     lhs_local_vars = [var for var in local_variables if var in write["line"].split("=")[1]]
                #     no_unsafe_writes= no_unsafe_writes and lhs_local_vars == []

                is_safe = no_unsafe_writes
                if not no_unsafe_writes:
                    if_nums = []
                    nest_nums = []
                    for write in var_writes:
                        for line in analyse_lines:
                            if line[3] == write["line_num"]:
                                if_nums.append(line[4])
                                nest_nums.append(line[1])
                    all_are_in_the_same_if = all([if_num == if_nums[0] for if_num in if_nums])
                    all_are_in_main_branch = all([nest_num == nest_nums[0] for nest_num in nest_nums])
                    is_safe = all_are_in_the_same_if or all_are_in_main_branch
                if is_safe:
                    safe_vars_to_inline.append(var)
            print("SAFE VARS TO INLINE",safe_vars_to_inline)
            variables = merge_dictionaries(self.static_variables, local_variables)
            possible_values = {}
            remove_indexes = []
            if_num = 0
            for line_index, line in enumerate(lines):
                line = line.strip()
                if "if" in line or "else" in line:
                    if_num += 1

                res_line = ""
                last_index = 0
                line = self.replace_segments(get_var_segments_in_line(line,variables),line,self.inline_replacer,local_variables,{
                    "possible_values": possible_values
                })
                lines[line_index] = line
                rhs_var = self.get_rhs_variable(line)
                if rhs_var in local_variables:
                    rhs_dim = len(variables[rhs_var]["dims"])
                    _ , num_of_looped_dims = get_dims_from_indexes(get_indexes(get_rhs_segment(line),rhs_var,rhs_dim),rhs_var)
                    if rhs_var in local_variables and rhs_dim <= 1 and rhs_dim == num_of_looped_dims and rhs_var in safe_vars_to_inline:
                        remove_indexes.append(line_index)
                        possible_values[rhs_var] = line.split("=")[1].replace("\n","")
            lines = [x[1] for x in enumerate(lines) if x[0] not in remove_indexes]
            lines = [line for line in lines if line != ""]
            return lines 

    def eliminate_while(self,lines):

        for line in lines:
            if "if (f /) then" in line:
                print("I fucked up")
                exit()
        for line in lines:
            if "if (f /) then" in line:
                print("I fucked up")
                exit()
        print("did first")
        #if there are some writes to flagged params then not safe to substitue
        writes = self.get_writes([(line,0) for line in lines])
        for flag_mapping in self.flag_mappings:
                if len([write for write in writes if write["variable"] == flag_mapping]) > 0:
                    del self.flag_mappings[flag_mapping]
                    break
        for line_index in range(len(lines)):
            for flag_mapping in self.flag_mappings:
                lines[line_index] = replace_variable(lines[line_index],flag_mapping,self.flag_mappings[flag_mapping])
                if "if (f /) then" in lines[line_index]:
                    print(flag_mapping,self.flag_mappings[flag_mapping])
                    print(line)
                    print(lines[line_index])
                    print("I fucked up mapping")
                    exit()
        # done twice since can have lflag_a -> lflag_b .or. lflag_c
        for line_index in range(len(lines)):
            for flag_mapping in self.flag_mappings:
                lines[line_index] = replace_variable(lines[line_index],flag_mapping,self.flag_mappings[flag_mapping])
                if "if (f /) then" in lines[line_index]:
                    print("I fucked up")
                    exit()


        print("done flag mapping")
        for line in lines:
            if "if (f /) then" in line:
                print("I fucked up")
                exit()
        print("passed flag map test")
        lines = self.transform_case(lines)
        print("passed transform case")
        local_variables = {parameter:v for parameter,v in self.get_variables([(line,0) for line in lines], {},self.file).items() }
        self.evaluate_ifs(lines,local_variables)
        print("passed evaluate ifs")
        writes = self.get_writes([(line,0) for line in lines])
        orig_lines = lines.copy()
        orig_lines.append("one more")
        print("start while")
        while(len(orig_lines) > len(lines)):
            print("start while")
            lines = self.eliminate_dead_branches(lines)
            for line in lines:
                if "if (f /) then" in line:
                    print("I fucked up")
                    exit()
            print("passed eliminate")
            writes = self.get_writes([(line,0) for line in lines])
            self.try_to_deduce_if_params(lines,writes,local_variables)
            self.evaluate_ifs(lines,local_variables)
            for line in lines:
                if "if (f /) then" in line:
                    print("I fucked up")
                    exit()
            print("passed eval")

            orig_lines = lines.copy()

            lines = self.eliminate_dead_branches(lines)
            for line in lines:
                if "if (f /) then" in line:
                    print("I fucked up")
                    exit()
            print("passed elim")
            writes = self.get_writes([(line,0) for line in lines])
            self.try_to_deduce_if_params(lines,writes,local_variables)
            self.evaluate_ifs(lines,local_variables)
            for line in lines:
                if "if (f /) then" in line:
                    print("I fucked up")
                    exit()
            print("passed eval")
        return lines
    def transform_line_stencil(self,line,num_of_looped_dims, local_variables, array_segments_indexes,rhs_var,vectors_to_replace):
        variables = merge_dictionaries(self.static_variables, local_variables)
        last_index = 0
        res_line = ""
        print("HMM",line)
        print("\n\n")
        print(array_segments_indexes)
        print("\n\n")
        #for testing remove later
        if "prof_lnt" in line:
            return ""
        make_vector_copies = False
        #whole nx writes become scalar writes 
        #3 dim writes are vectors
        if (
            rhs_var in local_variables and
                (
                    (num_of_looped_dims == 0 and len(local_variables[rhs_var]["dims"]) == 0 ) 
                    or (num_of_looped_dims == 1 and self.evaluate_indexes(local_variables[rhs_var]["dims"][0]) in ["nx","3"])
                ) 
                or (rhs_var in ["df","f"] or rhs_var in vectors_to_replace)
            ): 
            for i in range(len(array_segments_indexes)):
                segment = array_segments_indexes[i]
                if segment[0] in local_variables:
                    src = local_variables
                else:
                    src = self.static_variables
                if segment[0] == "df" or segment[0] == "f":
                    indexes = [self.evaluate_indexes(index) for index  in get_segment_indexes(segment,line, len(src[segment[0]]["dims"]))]
                    if not all([okay_stencil_index(x[1],x[0]) for x in enumerate(indexes[:-1])]):
                        print("How how to handle stencil indexes?")
                        print(line)
                        print(line[segment[1]:segment[2]])
                        print(indexes)
                        exit()
                    if ":" in indexes[-1]:
                        if is_vector_stencil_index(indexes[-1]):
                            make_vector_copies = True
                        else:
                            print("range in df index 3")
                            print(line[segment[1]:segment[2]])
                            print(indexes)
                            exit()
                    if segment[0] == "df":
                        vtxbuf_name = get_vtxbuf_name_from_index("DF_", indexes[-1])
                    elif segment[0] == "f":
                        vtxbuf_name = get_vtxbuf_name_from_index("VTXBUF_OUT_",indexes[-1])
                    res = vtxbuf_name
                else:
                    var_dims = src[segment[0]]["dims"]
                    indexes = [self.evaluate_indexes(index) for index  in get_indexes(line[segment[1]:segment[2]],segment[0],0)]
                    is_profile = segment[0] in self.static_variables and var_dims in [["nx"],["mx"],["ny"],["my"],["nz"],["mz"],["nx","my"]]
                    if is_profile:
                        #PROFILE_X
                        if var_dims in [["nx"],["mx"]]:
                            if indexes not in [[],["l1:l2"],["l1+3:l2+3"]]:
                                print(indexes)
                                print("WEIRD INDEX in profile_x")
                                print(indexes)
                                print(line[segment[1]:segment[2]])
                                exit()
                            profile_index = "0"
                            if indexes == ["l1+3:l2+3"]:
                                profile_index = "3"
                        #PROFILE_Y
                        elif var_dims in [["ny"],["my"]]:
                            if indexes not in  [["m"]]:
                                print("WEIRD INDEX in profile_y")
                                print(indexes)
                                print(line[segment[1]:segment[2]])
                                exit()
                            profile_index = "0"
                        #PROFILE_Z
                        elif var_dims in [["nz"],["mz"]]:
                            print("PROFILE_Z: ",segment[0])
                            if indexes not in [["n"],["n+3"]]:
                                print("WEIRD INDEX in profile_z")
                                print(indexes)
                                print(line[segment[1]:segment[2]])
                                exit()
                            profile_index = "0"
                            if indexes == ["n+3"]:
                                profile_index = "3"
                        #PROFILE_XY
                        elif var_dims in [["nx","my"]]:
                            if indexes not in [[":","m"]]:
                                print("WEIRD INDEX in profile_xy")
                                print(indexes)
                                print(line[segment[1]:segment[2]])
                                exit()
                            profile_index = "0"
                        else:
                            print("add profile mapping to profile", self.profile_mappings[segment[0]])
                            exit()
                        res = "AC_PROFILE_" + segment[0].upper() + "[" + profile_index + "]"
                    elif segment[0] in local_variables and self.evaluate_indexes(src[segment[0]]["dims"][0]) == "nx" and indexes in [[],[":"]]:
                        print("nx var -> scalar")
                        print(line[segment[1]:segment[2]])
                        print(segment[0])
                        print(line)
                        res = segment[0]
                    #global vec
                    elif segment[0] in self.static_variables and src[segment[0]]["dims"] == ["3"] and indexes in [["1"],["2"],["3"]]:
                        #do nothing
                        return line[segment[1]:segment[2]]
                    elif segment[0] in vectors_to_replace:
                        print("should make vec access")
                        res = line[segment[1]:segment[2]]
                        res = res.replace("(:,1)",".x")
                        res = res.replace("(:,2)",".y")
                        res = res.replace("(:,3)",".x")
                    #constant local array
                    elif len(src[segment[0]]["dims"]) == 1 and src[segment[0]]["dims"][0].isnumeric() and len(indexes) == 1:
                        return line[segment[1]:segment[2]]
                    else:
                        print("what to do?")
                        print(line[segment[1]:segment[2]])
                        print(segment[0])
                        print(line)
                        print("is static: ",segment[0] in self.static_variables)
                        print("is local: ",segment[0] in local_variables)
                        print("var dims",var_dims)
                        print(src[segment[0]])
                        print("make vector copies = ",make_vector_copies)
                        #for time being skip
                        if "p%uij" in segment[0]:
                            res = ""
                        exit()
                res_line = res_line + line[last_index:segment[1]]
                res_line = res_line + res 
                last_index = segment[2]
            res_line = res_line + line[last_index:]
            res_line = res_line + ";"
            if "%" in res_line:
                res_line = res_line.replace(";","")
                res_final_line = ""
                last_index = 0
                struct_segs = self.get_struct_segments_in_line(res_line,variables)
                print("\n\n")
                print("stuct segs",struct_segs)
                print("\n\n")
                for seg in struct_segs:
                    var_name,field = [part.strip() for part in seg[0].split("%",1)]
                    if var_name in local_variables:
                        src = local_variables
                    elif var_name in self.static_variables:
                        src = self.static_variables
                    if src[var_name]["type"] != "pencil_case":
                        print("what to do non pencil_case struct ?")
                        print("struct seg", seg[0], res_line[seg[1]:seg[2]])
                        exit()
                    else:
                        indexes = [self.evaluate_indexes(index) for index in get_segment_indexes(seg, res_line, 0)]
                        #replace 1:n -> with : if an index dim is n
                        for i,index in enumerate(indexes):
                            if ":" in index and index != ":":
                                lower,upper = [part.strip() for part in index.split(":")]
                                if lower == "1" and upper == self.struct_table[src[var_name]["type"]][field]["dims"][i]:
                                    indexes[i] = ":"
                        ##Pencil case will become a x dimensional profile.
                        ##If vec make three x dimensional profiles
                        if "(" not in res_line[seg[1]:seg[2]] or indexes == [":"] or (len(indexes) == 2 and indexes[0] == ":" and (indexes[1].isnumeric() or indexes[1] == ":")):
                            if make_vector_copies and (self.struct_table[src[var_name]["type"]][field]["dims"] == ["nx","3"]):
                                res = "AC_PROFILE_" + field.upper() + "$"
                            else:
                                res = "AC_PROFILE_" + field.upper()
                                if len(indexes) == 2:
                                    if indexes[1] == "1":
                                        res += "X"
                                    elif indexes[1] == "2":
                                        res += "X"
                                    elif indexes[1] == "3":
                                        res += "Y"
                                    else:
                                        print("WHAA")
                                        print(res_line[seg[1]:seg[2]])
                                        exit()
                        else:
                            print("weird array access in pencil case")
                            #for the time being skip
                            if "p%bij":
                                return ""
                            print(res_line[seg[1]:seg[2]])
                            print(indexes)
                            print(line)
                            exit()
                    res_final_line= res_final_line + res_line[last_index:seg[1]]
                    res_final_line = res_final_line + res 
                    last_index = seg[2]
                res_final_line= res_final_line + res_line[last_index:]
                print("res_final_line",res_final_line)
                res_line = res_final_line + ";"

            return res_line

        else:
            print("NO case for",line)
            print("RHS VAR",rhs_var)
            print(rhs_var in local_variables)
            print(local_variables[rhs_var]["dims"])
            exit()
    def transform_line_boundcond(self,line,num_of_looped_dims, local_variables, array_segments_indexes,rhs_var):
        last_index = 0
        res_line = ""
        if (num_of_looped_dims==0 or num_of_looped_dims==2) and (rhs_var in local_variables or rhs_var == "f"): 
            for i in range(len(array_segments_indexes)):
                segment = array_segments_indexes[i]
                print("seg",line[segment[1]:segment[2]])
                if segment[0] in local_variables:
                    src = local_variables
                else:
                    src = self.static_variables
                if segment[0] == "f":
                    indexes = [self.evaluate_indexes(index) for index in get_segment_indexes(segment,line, len(src[segment[0]]["dims"]))]
                    print("evaluated indexes",indexes)
                    vtxbuf_name = get_vtxbuf_name_from_index("VTXBUF_", indexes[-1])
                    new_var = f"vba.in[{vtxbuf_name}]"
                    indexes= indexes[:-1]
                    for i,index in enumerate(indexes):
                        if ":" in index:
                            indexes[i] = self.map_to_new_index(indexes[i],i,local_variables,line)
                        else:
                            #convert from 1 to 0 based index
                            indexes[i] = f"{indexes[i]}-1"
                    res = f"{new_var}[DEVICE_VTXBUF_IDX({','.join(indexes)})]"
                ##Value local to kernel i.e. from the viewpoint of a thread a scalar
                elif segment[0] in local_variables and len(local_variables[segment[0]]["dims"]) == 2 and num_of_looped_dims == 2:
                    res = segment[0] 
                else:

                    indexes = [self.evaluate_indexes(index) for index in get_segment_indexes(segment,line,len(src[segment[0]]["dims"]))]
                    num_of_looping_dim = -1
                    for i,index in enumerate(indexes):
                        if ":" in index:
                            num_of_looping_dim += 1
                        if index == ":":
                            possible_index = src[segment[0]]["dims"][i]
                            if possible_index == ":":
                                possible_index = local_variables[rhs_var]["dims"][num_of_looping_dim]
                            indexes[i] = self.map_to_new_index("1:" + possible_index,local_variables)
                        elif ":" in index:
                            indexes[i] = self.map_to_new_index(indexes[i],local_variables)
                        else:
                            indexes[i] = index
                    res = segment[0]
                    for index in indexes:
                        res = res + f"[{index}]"
                res_line = res_line + line[last_index:segment[1]]
                res_line = res_line + res 
                if ":" in res or "explicit_index" in res or ("(" in res and "DEVICE" not in res) or "'" in res:
                    print("NEED to transform more earlier")
                    print("LINE",line)
                    print("RES_LINE",res_line)
                    print("RES",res)
                    print("indexes",)
                    print("indexes",indexes)
                    print("seg",segment)
                    print("seg line",line[segment[1]:segment[2]])
                    exit()
                last_index = segment[2]
            res_line = res_line + line[last_index:]
            res_line = res_line + ";"
            if ":" in res_line or "explicit_index" in res_line or not has_balanced_parens(res_line) or ";" not in res_line or "()" in res_line:
                print("NEED to transform more")
                print("LINE",line)
                print("RES_LINE",res_line)
                print("RES",res)
                exit()
            return res_line

        else:
            print("NO case for",line)
            exit()
    def unroll_constant_loops(self,lines,local_variables):
        found_constant_loop = True
        while(found_constant_loop):
            lines,found_constant_loop = self.unroll_constant_loops_once(lines,local_variables)
        return lines
    def unroll_constant_loops_once(self,lines,local_variables):
        constant_loops_indexes = []
        in_constant_loop = False
        found_constant_loop = False
        loop_start_index  = 0
        for line_index, line in enumerate(lines):
            if not found_constant_loop:
                if in_constant_loop:
                    constant_loops_indexes.append(line_index)
                if line[:2] == "do" and "while" not in line and "end" not in line:
                    loop_index = self.get_writes_from_line(line,0)[0]["variable"]
                    lower,upper= [part.strip() for part in line.split("=")[1].split(",",1)]
                    if lower.isnumeric() and upper.isnumeric():
                        in_constant_loop = True
                        constant_loops_indexes = []
                        replacements = []
                        loop_start_index = line_index
                        add_replacement = int(lower)
                        while add_replacement <= int(upper):
                            replacements.append(str(add_replacement))
                            add_replacement += 1

                if "do" in line and "end" in line and in_constant_loop:
                    in_constant_loop = False
                    constant_loops_indexes.pop()
                    found_constant_loop = True
                    lines_to_add = []
                    remove_indexes = [loop_start_index,line_index]
                    # for replacement in replacements:
                    #     for x in lines_to_unroll:
                    #         lines_to_add.append(replace_variable(x,loop_index,replacement))
                    # for x in lines_to_add:
                    #     res_lines.append(x)
        if not found_constant_loop:
            return (lines,False)
        res_lines = []
        unrolled_lines = []
        for replacement in replacements:
            for x in constant_loops_indexes:
                unrolled_lines.append(replace_variable(lines[x],loop_index,replacement))
        for line_index,line in enumerate(lines):
            if line_index not in constant_loops_indexes and line_index not in remove_indexes:
                res_lines.append(line)
            if line_index == constant_loops_indexes[0]:
                res_lines.extend(unrolled_lines)
        print("res lines")
        for line in res_lines:
            print(line)
        return (res_lines,found_constant_loop)
    def transform_line(self,i,lines,local_variables,loop_indexes,symbol_table,initialization_lines,orig_params,transform_func,vectors_to_replace):
        line = lines[i]
        variables = merge_dictionaries(self.static_variables, local_variables)
        if is_init_line(line):
            vars_in_line = {}
            self.get_variables([(line,0)],vars_in_line, self.file)
            vars_to_declare = []
            for var in vars_in_line:
                if var.strip() in local_variables and (len(local_variables[var]["dims"]) == 0 or len(local_variables[var]["dims"]) == 2) and var.strip() not in orig_params:
                    vars_to_declare.append(var)
            if len(vars_to_declare) > 0 and local_variables[vars_to_declare[0]]["type"] != "pencil_case":
                return translate_to_c(local_variables[vars_to_declare[0]]["type"]) + " " + ", ".join(vars_to_declare) + ";"
            else:
                return ""
        if line.strip()  == "exit":
            return "continue"
        if "else" in line and "if" in line:
            return "}\n" + line.replace("then","{").replace("elseif","else if")
        if "else" in line:
            return "}\nelse {"
        if "if" in line and "then" in line:
            params = self.get_function_calls_in_line(line.replace("then",""),local_variables)[0]["parameters"]
            if len(params) == 1:
                return line.replace("then","{")
            else:
                print("what to do for if")
                print(line)
                exit()
        if "end" in line and "select" in line:
            return "}\n"
        if "select case" in line:
            select_case_var = line.split("(")[1].split(")")[0].lower()
            return f"switch({select_case_var})"+"{\n"
        if "case" in line and "default" in line:
            return "default:"
        if "case" in line:
            select_case_var = line.split("(")[1].split(")")[0]
            return f"case {select_case_var}:\n"
        if "end" in line and "do" in line:
            loop_indexes = loop_indexes[:-1]
            return "}\n"
        if "subroutine" in line and "end" in line:
            if self.test_to_c:
                return ""
            else:
                return "}\n"
        if "subroutine" in line:
            function_name = self.get_function_calls_in_line(line,local_variables)[0]["function_name"]
            if self.test_to_c:
                return ""
            else:
                return f"static __global__ void\n {function_name}(const int3 dims, VertexBufferArray vba)\n"+"{\n"
            
        if is_use_line(line):
            return ""
        if "do" in line[:2]:
            print("do line",line)
            loop_index = self.get_writes_from_line(line,0)[0]["variable"]
            lower,upper= [part.strip() for part in line.split("=")[1].split(",",1)]
            loop_indexes.append(loop_index)
            return f"for(int {loop_index} = {lower};{loop_index}<={upper};{loop_index}++)" +"{"
        if "endif" in line:
            return "}"
        if "endwhere" in line:
            return "}"
        original_line = line
        if new_is_variable_line(line):
            return ""

        func_calls = self.get_function_calls_in_line(line,local_variables)
        if len(func_calls) == 1 and func_calls[0]["function_name"] in ["der6","der"]:
            rest_params = func_calls[0]["parameters"][:2] + func_calls[0]["parameters"][3:]
            res = f"{func_calls[0]['parameters'][2]} = {func_calls[0]['function_name']}({','.join(rest_params)})"
            print("der call", line)
            print("res der call --->")
            print(res)
            return res
        if len(func_calls) > 0 and func_calls[0]["function_name"] in ["where"]:
            is_scalar_if = False
            print("where call",line)
            arr_segs_in_line = self.get_array_segments_in_line(line,variables)
            for seg in arr_segs_in_line:
                param_info = self.get_param_info((line[seg[1]:seg[2]],False),local_variables,self.static_variables)
                print("param info",param_info)
                is_scalar_if = is_scalar_if or (param_info[3] in [["nx","3"],["nx"]] )
            for seg in self.get_struct_segments_in_line(line,variables):
                param_info = self.get_param_info((line[seg[1]:seg[2]],False),local_variables,self.static_variables)
                print("param info",param_info)
                is_scalar_if = is_scalar_if or (param_info[3] in  [["nx"],["nx","3"]])
            if not is_scalar_if:
                print("what to about where")
                print(line)
                exit()
            print("return?")
            return line.replace("where","if",1) + "{"
        if "where" in line:
            print(func_calls)
            print("where still in line")
            print("where" in local_variables)
            print("where" in self.static_variables)
            print(line)
            exit()

        rhs_segment = get_rhs_segment(line)
        if rhs_segment is None:
            print("rhs seg is None")
            print("line: ",line)
            exit()
        rhs_var = self.get_rhs_variable(line)
        if rhs_var is None:
            print("rhs var is none")
            print(line)
            exit()
        rhs_var = rhs_var.lower()
        if rhs_var not in local_variables:
            print("WHAT TO DO rhs not in variables",line)
            print(rhs_var)

            #for the time being simply assume they are diagnostics writes so can simply remove them
            return ""
            print(rhs_var in self.static_variables)
            if rhs_var == "n":
                print("IS n_save in local_variables?","n_save" in local_variables)
            exit()
        dim = len(variables[rhs_var]["dims"])
        indexes = get_indexes(get_rhs_segment(line),rhs_var,dim)
        dims, num_of_looped_dims = get_dims_from_indexes(indexes,rhs_var)

        #line = self.transform_spread(line,[f":" for i in range(num_of_looped_dims)],local_variables)
        array_segments_indexes = self.get_array_segments_in_line(line,variables)
        return transform_func(line,num_of_looped_dims, local_variables, array_segments_indexes,rhs_var,vectors_to_replace)
    def transform_lines(self,lines,local_variables,transform_func):
        for i,line in enumerate(lines):
            if not has_balanced_parens(line):
                print("not balanced")
                print("--->")
                print(line)
                exit()
        for var in self.known_dims:
            self.static_variables[var]["dims"] = self.known_dims[var]
        symbol_table = {}
        #mapping that can be deduced for global flags
        lines = [line.replace("\n","") for line in lines]
        #transform cases to if statements to be able to analyse them together with other if cases

        orig_params = [param for param in self.get_function_calls_in_line(lines[0],local_variables)[0]["parameters"] if local_variables[param]["type"] != "pencil_case"]
        loop_indexes = []
        initialization_lines = []
        variables = merge_dictionaries(self.static_variables, local_variables)
        writes = self.get_writes([(line,0) for line in lines])
        lines = [line for line in lines if not any([func in line for func in self.safe_subs_to_remove])]
        ##Needed to remove size from variable dims
        lines = [self.expand_size_in_line(line,local_variables,writes) for line in lines]
        #get local variables back to get actual dims not size dims
        local_variables = {parameter:v for parameter,v in self.get_variables([(line,0) for line in lines], {},self.file).items() }
        #remove allocate and at the same time make sure the allocate is safe to remove i.e. refers to local variables
        remove_indexes = []
        for i,line in enumerate(lines):
            has_allocate = False 
            for func_call in self.get_function_calls_in_line(line,local_variables):
                if func_call["function_name"] == "allocate":
                    remove_indexes.append(i)
                    for param in [param.split("(")[0].strip() for param in func_call["parameters"] if "=" not in param and ("stat" not in param or "src" not in param or "source" not in param)]:
                        if param not in local_variables:
                            print("Subroutine allocates global variable", param)
                            print("So can't generate cuda code for it")
                            exit()

        lines = [x[1] for x in enumerate(lines) if x[0] not in remove_indexes]
        # lines = self.eliminate_while(lines)

        file = open(f"res-eliminated.txt","w")
        for line in lines:
            file.write(f"{line}\n")
        file.close()
        print("check eliminated file")


        global_loop_lines,iterators = self.get_global_loop_lines(lines,local_variables)
        if len(global_loop_lines) > 0:
            lines = self.remove_global_loops(lines,local_variables,global_loop_lines,iterators)

        #unroll forall lines
        for line_index,line in enumerate(lines):
            search_line = line.replace("forall","call forall")
            func_calls = self.get_function_calls_in_line(search_line,local_variables)
            if len(func_calls) == 1 and func_calls[0]["function_name"] == "forall" and len(func_calls[0]["parameters"]) == 1:
                write = self.get_writes_from_line((func_calls[0]["parameters"])[0],0)[0]
                iterator = write["variable"]
                if ":" not in write["value"]:
                    print("should unroll for all")
                    print("don't know how to do it")
                    print(line) 
                    exit()
                replacement_lower, replacement_upper = [part.strip() for part in write["value"].split(":")]

                res = self.replace_segments(self.get_array_segments_in_line(search_line,variables),search_line,self.global_loop_replacer,local_variables,{
                    "iterator": iterator, 
                    "replacement_lower": replacement_lower,
                    "replacement_upper": replacement_upper,
                    "replacement_index": 3,
                })
                save_var = False 
                #copypaste
                search_index = 0
                if "%" in res:
                    search_var = ""
                    save_var = False
                    buffer = ""
                    for i,char in enumerate(res):
                        print(char,buffer)
                        if char == "%":
                            print("TRUe",char)
                            save_var = True
                        if not(char.isalpha() or char.isnumeric()) and char != "%":
                            if save_var:
                                search_var = buffer
                                search_index = i
                                save_var = False
                            buffer = ""
                        else:
                            buffer = buffer + char
                    if save_var:
                        search_var = buffer
                        search_index = None
                    search_var = search_var.strip()
                    last_index = 0
                    res_final = ""
                    struct_segs = [seg for seg in get_variable_segments(res,search_var) if seg[0] != ""]
                    for seg in struct_segs:
                        var_name,field = [part.strip() for part in seg[0].split("%",1)]
                        if var_name in local_variables:
                            src = local_variables
                        elif var_name in self.static_variables:
                            src = self.static_variables
                        if src[var_name]["type"] != "pencil_case":
                            print("what to do non pencil_case struct ?")
                            print("struct seg", seg[0], res_line[seg[1]:seg[2]])
                            exit()
                        dims = self.struct_table[src[var_name]["type"]][field]["dims"]
                        indexes = get_segment_indexes(seg, res, 0)
                        for i,index in enumerate(indexes):
                            if iterator in index:
                                new_lower= index.replace(iterator,replacement_lower)
                                new_upper= index.replace(iterator,replacement_upper)
                                if new_lower == "1" and new_upper == dims[i]:
                                    indexes[i] = ":"
                                else:
                                    indexes[i] = new_lower + ":" + new_upper
                        res_final = res_final + res[last_index:seg[1]]
                        res_final= res_final + build_new_access(seg[0],indexes)
                        last_index = seg[2]
                    res_final = res_final + res[last_index:]
                    res = res_final

                func_call = self.get_function_calls_in_line(res,local_variables)[0]
                res = self.replace_func_call(res,func_call,"").replace("call","")
                lines[line_index] = res

        variables = merge_dictionaries(self.static_variables, local_variables)
        print("Unrolling constant loops")
        lines = [line.strip() for line in lines]
        #Unroll constant loops
        lines = self.unroll_constant_loops(lines,local_variables)
        file = open("res-unroll.txt","w")
        for line in lines:
            file.write(f"{line}\n")
        file.close()
        print("Check unroll file")

        print("replace 1d vecs with 3 1d arrays")
        vectors_to_try_to_replace = []
        for var in local_variables:
            vectors_to_try_to_replace.append(var)
        print("dline_1 is in self.static","dline_1" in self.static_variables)
        #Make sure are actually vectors
        vectors_to_replace = []
        for vector in vectors_to_try_to_replace:
            if self.is_vector(lines, vector, local_variables):
                vectors_to_replace.append(vector)
        vectors_to_replace.append("dline_1")

        lines = [line.strip() for line in lines if line.strip() != ""]
        file = open("res-inlined.txt","w")
        for line in lines:
            file.write(f"{line}\n")
        file.close()
        print("Check inlined file")
        #transform spreads into do loops
        res_lines = []
        for line_index, line in enumerate(lines):
            spread_calls= [call for call in self.get_function_calls_in_line(line,variables) if call["function_name"] == "spread"]
            if len(spread_calls) > 1:
                print("multiple spread calls",line)
                exit()
            elif len(spread_calls) == 1:
                call = spread_calls[0]
                print("spread call")
                print(line)
                lhs = call["parameters"][0]
                redundant_index = call["parameters"][1]
                rhs_var = self.get_rhs_variable(line)
                rhs_info = self.get_param_info((rhs_var,False),local_variables,self.static_variables)
                res_line = self.replace_func_call(line,call,lhs)
                rhs_segment = self.get_array_segments_in_line(line,variables)[0]
                print(rhs_segment)
                rhs_indexes = get_segment_indexes(rhs_segment,res_line,0)
                if rhs_indexes == [] and len(variables[rhs_var]["dims"]) == 1:
                    new_rhs = f"{rhs_var}(:)"
                    res_line = new_rhs + res_line[rhs_segment[2]:]
                    res_lines.append(f"do explicit_index = 1,{rhs_info[3][0]}")
                    res_lines.append(res_line)
                    res_lines.append("enddo")
                else:
                    print("have to append it to")
                    exit()
            
            else:
                res_lines.append(line)
        lines = res_lines
        #try to transform maxvals into max before transform 
        for line_index,line in enumerate(lines):
            maxval_calls = [call for call in self.get_function_calls_in_line(line,variables) if call["function_name"] == "maxval"]
            if len(maxval_calls) > 1:
                print("multiple maxval calls")
                print(line)
                exit()
            elif len(maxval_calls) == 1:
                print("maxvall call",maxval_calls[0])
                first_param_info = self.get_param_info((maxval_calls[0]["parameters"][0],False),local_variables,self.static_variables)
                print("first param info", first_param_info)
                #max of nx (possibly vector components), is safe to take max
                if ((len(maxval_calls[0]["parameters"]) == 1) or (len(maxval_calls[0]["parameters"]) == 2 and maxval_calls[0]["parameters"][1] in ["dim=2","2"])) and first_param_info[3] in [["nx"], ["nx","3"]]:
                    print("maxval before",line)
                    lines[line_index] = self.replace_func_call(line,maxval_calls[0],f"max({maxval_calls[0]['parameters'][0]})")
                    print("maxval after",lines[line_index])
                else:
                    print(line)
                    print("first param info",first_param_info)
                    print("no case maxval")
                    exit()
        #replace all calls to tiny with predefined float 
        for line_index,line in enumerate(lines):
            for call in [x for x in self.get_function_calls_in_line(line,variables) if x["function_name"] == "tiny"]:
                type = self.get_param_info((call["parameters"][0],False), local_variables,self.static_variables)[2]
                if type == "real":
                    line = self.replace_func_call(line,call,"AC_tiny_val")
                    lines[line_index] = line
                else:
                    print("what to do with this tiny type?",type)
                    exit()
        #alog -> log
        for line_index,line in enumerate(lines):
                alog_calls =  [x for x in self.get_function_calls_in_line(line,variables) if x["function_name"] == "alog"]
                while(len(alog_calls)>0):
                    call = alog_calls[0]
                    line = self.replace_func_call(line,call,f"log({','.join(call['parameters'])})")
                    lines[line_index] = line
                    alog_calls =  [x for x in self.get_function_calls_in_line(line,variables) if x["function_name"] == "alog"]
        #for calls [real] just return the param
        for line_index,line in enumerate(lines):
                real_calls = [x for x in self.get_function_calls_in_line(line,variables) if x["function_name"] == "real"]
                while(len(real_calls)>0):
                    call = real_calls[0]
                    if len(call["parameters"]) == 1:
                        print("real before",line)
                        line = self.replace_func_call(line,call,call['parameters'][0])
                        print("real after",line)
                        lines[line_index] = line
                        real_calls = [x for x in self.get_function_calls_in_line(line,variables) if x["function_name"] == "real"]
                    else:
                        print("multiple params in real")
                        print(line)
                        exit()

        #transform sum calls
        writes = self.get_writes([(line,0) for line in lines])
        self.try_to_deduce_if_params(lines,writes,local_variables)
        for i,line in enumerate(lines):
            res = self.transform_line(i,lines,local_variables,loop_indexes,symbol_table,initialization_lines,orig_params, transform_func,vectors_to_replace)
            lines[i] = res
        file = open("res-transform.txt","w")
        for line in lines:
            file.write(f"{line}\n")
        file.close()
        lines = [line.replace(".false.","false") for line in lines]
        lines = [line.replace(".true.","true") for line in lines]
        lines = [line.replace("/=","!=") for line in lines]

        lines = [replace_exp(line) for line in lines]
        file = open("res-replace_exp.txt","w")
        for line in lines:
            file.write(f"{line}\n")
        file.close()
        print("Check replace exp file")
 
        #for testing vtxbuf_ss -> vtxbuf_entropy
        lines = [line.replace("VTXBUF_SS","VTXBUF_ENTROPY") for line in lines]
        #close empty defaults
        for i,line in enumerate(lines):
            if "default:" in line and i<len(lines)-1 and "}" in lines[i+1]:
                lines[i] = line + ";"
        print("Check transform file")
        res_lines = []
        for line in lines:
            res_lines.extend([x for x in line.split("\n") if x != ""])
        lines = res_lines
        static_vars = []
        for i,line in enumerate(lines):
            static_variables_in_line= list(set([var.lower() for var in get_used_variables_from_line(line) if var.lower() in self.static_variables]))
            static_vars.extend(static_variables_in_line)
        static_vars = list(set(static_vars))
        file = open("ac_declarations.h","w")
        for var in static_vars:
            type = ""
            if self.static_variables[var]["type"] == "real":
                type = "real"
            if self.static_variables[var]["type"] == "integer":
                type = "int"
            if self.static_variables[var]["type"] == "logical":
                type = "int"
            if len(self.static_variables[var]["dims"]) == 1:
                type = type + "Array"
            if len(self.static_variables[var]["dims"]) == 2:
                type = type + "2dArray"
            if "intArray" not in type and type != "" and type != "Array":
                file.write(f"{type} AC_{var}\n")
        file.close()



        lines = self.replace_var_in_lines(lines, "nx","DCONST(AC_nx).x")
        lines = self.replace_var_in_lines(lines, "ny","DCONST(AC_ny).y")
        lines = self.replace_var_in_lines(lines, "nz","DCONST(AC_nz).z")

        lines = self.replace_var_in_lines(lines, "mx","(DCONST(AC_nx)+2*NGHOST)")
        lines = self.replace_var_in_lines(lines, "my","(DCONST(AC_ny)+2*NGHOST)")
        lines = self.replace_var_in_lines(lines, "mz","(DCONST(AC_nz)+2*NGHOST)")

        lines = self.replace_var_in_lines(lines, "l1","DCONST(AC_nx_min)")
        lines = self.replace_var_in_lines(lines, "l2","DCONST(AC_nx_max)")

        lines = self.replace_var_in_lines(lines, "m1","DCONST(AC_ny_min)")
        lines = self.replace_var_in_lines(lines, "m2","DCONST(AC_ny_max)")

        lines = self.replace_var_in_lines(lines, "n1","DCONST(AC_nz_min)")
        lines = self.replace_var_in_lines(lines, "n2","DCONST(AC_nz_max)")

        orig_lines = lines.copy()
        for i,line in enumerate(lines):
            static_variables_in_line= list(set([var for var in get_used_variables_from_line(line) if var.lower() in self.static_variables]))
            res_line = line
            for var in static_variables_in_line:
                ##all uppercase means that it is a profile
                #pow is a function in this context not a parameter
                if var.lower() not in ["bot","top","nghost","pow"] and var.lower() not in local_variables and var.upper() != var:
                    res_line = replace_variable(res_line, var, f"DCONST(AC_{var.lower()})")
            lines[i] = res_line

        if self.test_to_c:
            idx_line = "const int3 idx = {i, j, k};"
            lines = [line.replace("vba.in","mesh_test.vertex_buffer") for line in lines]
            res_lines = [idx_line] + lines
            print("TEST TO C")
            return res_lines
        static_vars = []
        rest_params = ""
        for param in orig_params:
            print("in orig params",param)
            rest_params = rest_params + "," + translate_to_c(local_variables[param]["type"]) + " " + param
        lines[1] = lines[1].replace("[rest_params_here]",rest_params)
        lines = [line.replace("nghost","NGHOST") for line in lines]
        file = open("res-3.txt","w")
        for line in lines:
            file.write(f"{line}\n")
        file.close()
        for i,line in enumerate(lines):
            func_calls = self.get_function_calls_in_line(line,local_variables)
            if i>1:
                for func_call in func_calls:
                    if func_call["function_name"].lower() not in ["sqrt","abs","tanh","min","max","der","der6","pow","DEVICE_VTXBUF_IDX".lower(),"DCONST".lower(),"exp","log","if","else","for","sin","cos"] and not "[" in func_call["function_name"]:
                        print("STILL FUNC CALLS in line:",line,i)
                        print(func_call)
                        exit()
        vertexIdx_line = "const int3 vertexIdx = (int3){\nthreadIdx.x + blockIdx.x * blockDim.x,\nthreadIdx.y + blockIdx.y * blockDim.y,\nthreadIdx.z + blockIdx.z * blockDim.z,\n};\n"
        check_line = "if (vertexIdx.x >= dims.x || vertexIdx.y >= dims.y || vertexIdx.z >= dims.z) {\nreturn;\n}\n"
        idx_line = "const int3 idx = {vertexIdx.x, vertexIdx.y, vertexIdx.z};"

        declarations_line = ""
        for var in list(set(initialization_lines)):
            print("in dec lines",var)
            declarations_line = declarations_line + translate_to_c(local_variables[var.lower()]["type"]) + " " + var + ";\n"
        res_lines = lines[:3] + [declarations_line, vertexIdx_line, check_line, idx_line] +lines[3:]
        lines = res_lines
        return lines
    def deduce_value(self,variable,writes,local_variables):
        # print("Trying to deduce value for", variable)
        var_writes = [write for write in writes if write["variable"] == variable and write["line"].split(" ")[0].strip() != "do"]
        if len(var_writes) > 1:
            if all([write["line"].split("=")[1].strip() == var_writes[0]["line"].split("=")[1].strip() for write in var_writes]):
                var_write = var_writes[0]
                value = var_write["line"].split("=")[1].replace("\n","").strip()
                local_variables[variable]["value"] = value
                self.known_values[variable] = value
                # print("deduced value for variable ",variable," = ",value)
                return
            return
        elif len(var_writes) == 1:
            var_write = var_writes[0]
            value = var_write["line"].split("=")[1].replace("\n","").strip()
            local_variables[variable]["value"] = value
            self.known_values[variable] = value
            # print("deduced value for variable ",variable," = ",value)
            return
        return
    
    def expand_size_in_line(self, line,local_variables,writes):
        func_calls = self.get_function_calls_in_line(line,local_variables,False)
        need_to_transform_size = any([func_call["function_name"] == "size" for func_call in func_calls])
        while need_to_transform_size:
            func_call = list(filter(lambda func_call: func_call["function_name"] == "size", func_calls))[0]
            replacement = self.get_size(func_call,local_variables,local_variables,writes)
            line = self.replace_func_call(line, func_call,replacement)
            func_calls =self.get_function_calls_in_line(line,local_variables,False)
            need_to_transform_size = any([func_call["function_name"] == "size" for func_call in func_calls])
        return line
    def get_global_loop_lines(self,lines,local_variables):
        global_loop_lines = []
        in_global_loop = False
        number_of_dos = 0
        number_of_enddos= 0
        lines_in_loop = []
        iterators = []
        for i,line in enumerate(lines):
            if in_global_loop:
                lines_in_loop.append((line,i))
            if "end" in line and "do" in line and in_global_loop:
                number_of_enddos += 1
                in_global_loop = not number_of_enddos == number_of_dos
                if not in_global_loop:
                    global_loop_lines.append(lines_in_loop)
                    lines_in_loop = []
                    number_of_dos = 0
                    number_of_enddos = 0
            elif self.offload_type == "stencil" and "do" in line and "1,nx" in line:
                write = self.get_writes_from_line(line,local_variables)[0]
                iterators.append((write["variable"],"1","nx",0))
                if in_global_loop:
                    print("Can't handle nested global loops")
                    exit()
                in_global_loop = True
                number_of_dos = 1
                number_of_enddos = 0
                lines_in_loop = []
                lines_in_loop.append((line,i))
            elif self.offload_type == "boundcond" and "do" in line and "1,my" in line:
                write = self.get_writes_from_line(line,local_variables)[0]
                iterators.append((write["variable"],"1","my",1))
                if in_global_loop:
                    print("Can't handle nested global loops")
                    exit()
                in_global_loop = True
                number_of_dos = 1
                number_of_enddos = 0
                lines_in_loop = []
                lines_in_loop.append((line,i))
        return (global_loop_lines,iterators)
    def try_to_deduce_if_params(self,lines,writes,local_variables):
        for var in local_variables:
            self.deduce_value(var,writes,local_variables)
        # for line in lines:
        #     check_line = line.replace("then","").replace("elseif","if").replace("if","call if")
        #     if_func_calls = list(filter(lambda func_call: func_call["function_name"] == "if", self.get_function_calls_in_line(check_line,local_variables)))
        #     if len(if_func_calls) == 1:
        #         if len(if_func_calls[0]["parameters"]) == 1:
        #             param = if_func_calls[0]["parameters"][0]
        #             params = [part.replace(".not.","").strip() for part in param.split(".and.")]
        #             for par in params:
        #                 for var in [var for var in local_variables if "value" not in local_variables[var]]:
        #                     if var == par:
        #                         self.deduce_value(var,writes,local_variables)

    def choose_correct_if_lines(self,lines):
        if_lines = []
        for line_index,line in enumerate(lines):
            if "if" in line or "else" in line:
                if_lines.append((line,line_index))
        possibilities = []
        found_true = False 
        for x_index, x in enumerate(if_lines):
            line = x[0]
            if "if" in line and ".false." in line and ".true." not in line and ".and." not in line and ".or." not in line:
                pass
            elif "if" in line and "else" not in line and ".true." in line and ".false." not in line and ".and." not in line and ".or." not in line:
                found_true = True
                possibilities= [x_index]
            elif not found_true and "end" not in line: 
                possibilities.append(x_index)
        print("possibilities",possibilities)
        if len(possibilities) == 1 and len(if_lines) > 2:
            correct_index = possibilities[0] 
            lower = if_lines[correct_index][1] + 1
            upper = if_lines[correct_index+1][1]
            return lines[lower:upper]
        elif len(possibilities) == 1:
            correct_index = possibilities[0] 
            line = lines[if_lines[correct_index][1]]
            if "if" in line and "else" not in line and "then" in line and ".true." in line and ".false." not in line and ".and." not in line and ".or." not in line:
                lower = if_lines[correct_index][1] + 1
                upper = if_lines[correct_index+1][1]
                return lines[lower:upper]
            elif "if" in line and "else" not in line and "then" in line and ".false." in line and ".true." not in line and ".and." not in line and ".or." not in line:
                return []
            else:
                return lines
        elif len(possibilities) == 0:
            return []
        return lines
        
    def eliminate_dead_branches(self,lines):
        orig_lines = lines.copy()
        orig_lines.append("one more")
        while(len(orig_lines) > len(lines)):
            orig_lines = lines.copy()
            print("\nStarting elimination\n")
            file = open("res-eliminated.txt","w")
            for line in lines:
                file.write(f"{line}\n") 
            file.close()
            lines = self.eliminate_dead_branches_once(lines)
            file = open("res-eliminated.txt","w")
            for line in lines:
                file.write(f"{line}\n") 
            file.close()
        return lines
    def get_analyse_lines(self,lines):
        if_num = 0
        analyse_lines = []
        remove_indexes = []
        done = False
        max_if_num = 0
        nest_num = 0
        case_num = 0
        if_nums = [0]
        for line_index,line in enumerate(lines):
            if "if" in line and "then" in line and "elif" not in line and "else" not in line:
                if_num = max_if_num + 1
                if_nums.append(if_num)
                max_if_num = max(if_num, if_num)
                nest_num += 1
            if "if" in line or "else" in line:
                case_num += 1
            analyse_lines.append((line,nest_num, if_num,line_index,case_num))
            if "if" in line and "end" in line:
                nest_num -= 1
                if_nums.pop()
                if_num = if_nums[-1]
        return analyse_lines
    def has_no_ifs(self,lines):
        return all([x[1] == 0 for x in self.get_analyse_lines(lines)])
    def eliminate_dead_branches_once(self,lines):
        if_num = 0
        analyse_lines = []
        remove_indexes = []
        done = False
        max_if_num = 0
        nest_num = 0
        if_nums = [0]

        for line_index,line in enumerate(lines):
            if "if" in line and "then" in line and "elif" not in line and "else" not in line:
                print("start",line)
                if_num = max_if_num + 1
                if_nums.append(if_num)
                max_if_num = max(if_num, if_num)
                nest_num += 1
            analyse_lines.append((line,nest_num, if_num,line_index))
            if "endif" in line or "end if" in line:
                print("end",line)
                nest_num -= 1
                if_nums.pop()
                if_num = if_nums[-1]

        for if_num in range(1,max_if_num+1):
            if_lines = [line for line in analyse_lines if line[2] == if_num and line[1] > 0]
            choices = []
            for line in if_lines:
                if ("if" in line[0] and ("then" in line[0] or "end" in line[0])) or "else" in line[0]:
                    choices.append(line)
            print("CHOICES")
            print(choices)
            possibilities = []
            found_true = False
            res_index = 0
            for choice_index, choice in enumerate(choices):
                line = choice[0]
                if "if" in line and ".false." in line and ".true." not in line and ".and." not in line and ".or." not in line:
                    pass
                elif "if" in line and ".true." in line and ".false." not in line and ".and." not in line and ".or." not in line:
                    found_true = True
                    possibilities = [choice]
                    res_index = choice_index
                elif not found_true and "end" not in line:
                    possibilities.append(choice)
                    res_index = choice_index
            print("POSSIBILITES")
            print(possibilities)
            if len(possibilities) == 0:
                starting_index = choices[0][3]
                ending_index = choices[-1][3]
                for index in range(starting_index,ending_index+1):
                    remove_indexes.append(index)
            if len(possibilities) == 1 and (len(choices) > 2 or found_true):
                starting_index = choices[0][3]
                ending_index = choices[-1][3]

                keep_starting_index = possibilities[0][3]+1
                #till next conditional or end
                keep_ending_index = choices[res_index+1][3]-1
                for index in range(starting_index, ending_index+1):
                    if not (index >= keep_starting_index and index<=keep_ending_index):
                        remove_indexes.append(index)
        new_lines = [x[1] for x in enumerate(lines) if x[0] not in remove_indexes]


        #choose correct single line ifs
        remove_indexes = []
        for line_index, line in enumerate(new_lines):
            if "if" in line and "else" not in line and "then" not in line and ".false." in line and ".true." not in line and ".and." not in line and ".or." not in line:
                remove_indexes.append(line_index)
            elif "if" in line and "else" not in line and "then" not in line and ".true." in line and ".else." not in line and ".and." not in line and ".or." not in line:
                new_lines[line_index] = new_lines[line_index].replace("if (.true.)","").replace("if(.true.)","")
        return [x[1] for x in enumerate(new_lines) if x[0] not in remove_indexes]
    
    def split_line(self, line):
        if line[0] == "!":
            return [line]
        lines = []
        split_indexes = []
        num_of_single_quotes = 0
        num_of_double_quotes = 0
        num_of_left_brackets = 0
        num_of_right_brackets = 0
        for iter_index in range(len(line)):
            if line[iter_index] == "'":
                num_of_single_quotes += 1
            if line[iter_index] == '"':
                num_of_double_quotes += 1
            if line[iter_index] == "(":
                num_of_left_brackets += 1
            if line[iter_index] == ")":
                num_of_right_brackets += 1
            if line[iter_index] == ";" and num_of_single_quotes%2 == 0 and num_of_double_quotes%2 == 0 and num_of_left_brackets == num_of_right_brackets:
                split_indexes.append(iter_index)
        start_index = 0
        for split_index in split_indexes:
            lines.append(line[start_index:split_index])
            start_index=split_index+1
        lines.append(line[start_index:])
        return filter(lambda x: x != "",lines)

    def get_function_declarations(self, lines, function_declarations, filename):
        ## Somewhat buggy, also returns variables in case of a header file but this is not a problem for the use case
        in_type = False
        for (line, count) in lines:
            parts = line.split("::")
            is_variable = True
            start =  parts[0].strip()
            type = start.split(",")[0].strip()
            if type == "public":
                is_variable = False
            if line.split(" ")[0] == "type" and len(line.split("::")) == 1:
                in_type = True
            if line.split(" ")[0] == "endtype":
                in_type = False
            if not is_variable and not in_type and len(parts)>1:
                allocatable = "allocatable" in [x.strip() for x in start.split(",")]
                saved_variable = "save" in [x.strip() for x in start.split(",")]
                public = "public" in [x.strip() for x in start.split(",")]
                dimension = re.search("dimension\s*\((.+?)\)", start)
                if dimension is None:
                    dimension = []
                else:
                    dimension = [x.strip() for x in dimension.group(1).split(",")]
                line_variables = parts[1].strip()
                variable_names = self.get_variable_names_from_line(line_variables)
                for i, variable_name in enumerate(variable_names):
                    function_declarations.append(variable_name)
        return function_declarations
    def expand_function_call_main(self, subroutine_lines, subroutine_name, filename, function_to_expand,global_init_lines,sub_modules):
        local_variables = {parameter:v for parameter,v in self.get_variables(subroutine_lines, {},filename).items() }
        function_calls = self.get_function_calls(subroutine_lines, local_variables)
        replaced_function_call_lines,is_function,res_type = self.expand_function_call(subroutine_lines, subroutine_name, filename, function_to_expand,local_variables,global_init_lines)
        subroutine_lines = [line for line in subroutine_lines if not is_init_line(line[0])]
        subroutine_lines = [line for line in subroutine_lines if not is_use_line(line[0])]
        init_var_names = []
        global_init_lines = list(set(global_init_lines))
        # for line in global_init_lines:
        #     new_init_var_names =  {parameter:v for parameter,v in self.get_variables([(line,0)], {},filename).items() }.keys()
        #     for new_var in new_init_var_names:
        #         if new_var in init_var_names:
        #             print("conflicting duplicate names",new_var)
        #             for line in global_init_lines:
        #                 if new_var in line:
        #                     print(line)
        #             exit()
        #         else:
        #             init_var_names.append(new_var)
        if is_function:
            if res_type == "":
                print("what went wrong?")
                print(function_to_expand)
                exit()
            subroutine_lines = subroutine_lines[:1] +  [(f"use {module}",0) for module in sub_modules] + [(line,0) for line in global_init_lines] + [(f"{res_type} :: {function_to_expand}_return_value",0)] + subroutine_lines[1:]
        else:
            subroutine_lines = subroutine_lines[:1] +  [(f"use {module}",0) for module in sub_modules] + [(line,0) for line in global_init_lines] + subroutine_lines[1:]
        has_replaced_call = False
        line_num = 0
        while(not has_replaced_call and line_num < len(subroutine_lines)-1):
            if f"{function_to_expand} " in subroutine_lines[line_num][0] or f"{function_to_expand}(" in subroutine_lines[line_num][0]:
                if is_function:
                    func_calls = self.get_function_calls_in_line(subroutine_lines[line_num][0],local_variables)
                    func_call = [func_call for func_call in func_calls if func_call["function_name"] == function_to_expand][0]
                    subroutine_lines[line_num] = (self.replace_func_call(subroutine_lines[line_num][0], func_call, f"{function_to_expand}_return_value"), subroutine_lines[line_num][1])
                    res_lines = subroutine_lines[:line_num] + replaced_function_call_lines +subroutine_lines[line_num:]
                else:
                    func_calls = self.get_function_calls_in_line(subroutine_lines[line_num][0],local_variables)
                    func_call = [func_call for func_call in func_calls if func_call["function_name"] == function_to_expand][0]
                    before_line = self.replace_func_call(func_call["line"],func_call,"$$").split("$$",1)[0].replace("call","").strip()
                    after_line = self.replace_func_call(func_call["line"],func_call,"$$").split("$$",1)[1].strip()
                    if "if" in before_line and after_line == "":
                        print("should add if to before_line")
                        print(func_call["line"])
                        before_line = before_line + " then"
                        after_line = "endif"
                    res_lines = subroutine_lines[:line_num] + [(before_line,0)] + replaced_function_call_lines + [(after_line,0)] + subroutine_lines[line_num+1:]
                subroutine_lines = res_lines
                has_replaced_call = True
            line_num = line_num + 1
        # new_function_calls = self.get_function_calls(subroutine_lines, local_variables)
        # for function_call in new_function_calls:
        #     if function_call["function_name"] != subroutine_name and function_call["function_name"] not in self.ignored_subroutines:
        #         if function_call["function_name"] == function_to_expand:
        #             print("DIDN'T i already expand this?",function_to_expand,function_call["parameters"])
        #             exit()
        #         subroutine_lines = self.expand_function_call_main(subroutine_lines, subroutine_name, filename, function_call["function_name"],depth+1)
        return subroutine_lines
    def populate_pencil_case(self):
        self.struct_table["pencil_case"] = {}
        for file in self.used_files:
            for (line,count) in self.get_lines(file,include_comments=True):
                if line[0] == "!" and "PENCILS PROVIDED" in line:
                    line = line.replace("!","").replace("PENCILS PROVIDED","")
                    fields = []
                    in_type = False
                    field = ""
                    for char in line:
                        if char == " " or char == "," or char == ";":
                            if not in_type:
                                fields.append(field)
                                field = ""
                            else:
                                field += char
                        elif char == "(":
                            in_type = True
                            field += char
                        elif char == ")":
                            in_type = False
                            fields.append(field)
                            field = ""
                        else:
                            field += char
                    fields.append(field)
                    fields = [field.strip().lower() for field in fields if field != ""]
                    for field in fields:
                        field_name = field.split("(")[0].strip()
                        field_dims = ["nx"]
                        if "(" in field:
                            field_dims += (field.split("(",1)[1].split(")")[0].split(","))
                        if field_name in self.struct_table["pencil_case"]:
                            if not len(self.struct_table["pencil_case"][field_name]["dims"]) == len(field_dims):
                                print("disagreeing pencils provided")
                                print("file", file)
                                print(field,self.struct_table["pencil_case"][field_name])
                                print(line)
                                exit()
                        else:
                            self.struct_table["pencil_case"][field_name] = {"type": "real", "dims": field_dims}

    def expand_function_call(self,lines,subroutine_name, filename, function_to_expand,variables_in_scope,global_init_lines):

        print(f"Expanding {function_to_expand} in {subroutine_name} in file {filename} inline_num {self.inline_num}")
        mpi_calls = ["mpi_send","mpi_barrier","mpi_finalize"]
        if function_to_expand in mpi_calls:
            print("MPI call not safe :(")
            exit()
        file_paths = self.find_subroutine_files(function_to_expand)
        #if file_paths is [] then the function is only present in the current file and not public
        if file_paths == []:
            file_paths = [filename]
        modules = self.get_used_modules(lines)
        self.parse_file_for_static_variables(filename)
        for module in self.get_subroutine_modules(filename, subroutine_name):
            self.parse_module(module)
        local_variables = {parameter:v for parameter,v in self.get_variables(lines, {},filename).items() }
        function_calls = self.get_function_calls(lines, local_variables)
        for function_call in function_calls:
            if function_call["function_name"] == function_to_expand:
                parameter_list = self.get_static_passed_parameters(function_call["parameters"],local_variables,self.static_variables)
                if len(file_paths)>1:
                    function_to_expand_filename = self.choose_right_module(file_paths)
                else:
                    function_to_expand_filename = file_paths[0]
                #return on the first call
                return self.get_replaced_body(function_to_expand, function_to_expand_filename,parameter_list, function_call,variables_in_scope,global_init_lines)
        print("HERE")
        return ([],False,None)
    def transform_case(self,lines):
        found_case = False 
        remove_indexes = []
        case_indexes = []
        case_param = ""
        case_params = []
        end_index = 0
        in_case = False
        for i, line in enumerate(lines):
            func_calls = self.get_function_calls_in_line(line,{})
            if in_case and len(func_calls) == 1 and func_calls[0]["function_name"] == "case":
                case_indexes.append(i)
                case_params.append(func_calls[0]["parameters"])
            if in_case and "default" in line and "case" in line:
                case_indexes.append(i)
            if "case" in line and "select" in line:
                in_case = True
                found_case = True
                remove_indexes = [i]
                case_indexes = []
                case_param = func_calls[0]["parameters"][0]
                case_params = []
            if "endselect" in line or "end select" in line:
                end_index = i
                in_case = False
                break
        if not found_case:
            return lines
        res_lines = []
        for i,x in enumerate(case_params):
            if i == 0:
                inside = ".or.".join([f"{case_param} == {y}" for y in x])
                res = f"if({inside}) then" 
            else:
                inside = ".or.".join([f"{case_param} == {y}" for y in x])
                res = f"else if({inside}) then" 
            if "if (f /) then" in res:
                print("I fucked up in transform case")
                exit()
            res_lines.append(res)
        #default case is handled separately
        res_lines.append("else")
        for j,i in enumerate(case_indexes):
            lines[i] = res_lines[j]
        lines[end_index] = "endif"
        lines = [x[1] for x in enumerate(lines) if x[0] not in remove_indexes]
        return self.transform_case(lines)
    def inline_all_function_calls(self,filename,subroutine_name,add_init_lines=False):
        self.parse_module(self.get_own_module(filename))
        self.load_static_variables(filename)
        new_lines = [(line[0].lower(),line[1]) for line in self.get_subroutine_lines(subroutine_name,filename)]
        
        init_lines = [line[0] for line in new_lines if is_init_line(line[0])]
        global_init_lines = init_lines
        global_init_lines = list(set(global_init_lines))
        elim_lines = self.eliminate_while([line[0] for line in new_lines])
        # elim_lines = [line[0] for line in new_lines]
        for line in elim_lines:
            if "if(f /) then" in line:
                print("I fucked up")
                exit()

        new_lines = [(line,0) for line in elim_lines]
        local_variables = {parameter:v for parameter,v in self.get_variables(new_lines, {},filename).items() }
        func_calls_to_replace = [call for call in self.get_function_calls(new_lines,local_variables) if call["function_name"] != subroutine_name and call["function_name"] not in self.ignored_subroutines]
        sub_modules = self.get_subroutine_modules(filename,subroutine_name)
        while len(func_calls_to_replace) != 0:
            new_lines = self.expand_function_call_main(new_lines, subroutine_name, filename,func_calls_to_replace[0]["function_name"],global_init_lines,sub_modules)
            local_variables = {parameter:v for parameter,v in self.get_variables(new_lines, {},filename).items() }
            modules =  [module_name for module_name in  self.get_used_modules(new_lines) if module_name not in self.ignored_modules]
            for module in modules:
                self.parse_module(module)
            func_calls_to_replace = [call for call in self.get_function_calls(new_lines,local_variables) if call["function_name"] != subroutine_name and call["function_name"] and call["function_name"] not in self.ignored_subroutines]
            file = open(f"res.txt","w")
            print("num_of_inlined", self.inline_num)
            for line in new_lines:
                file.write(f"{line[0]}\n")
            file.close()
            self.inline_num += 1
        if "inlined_lines" not in self.func_info[subroutine_name]:
            self.func_info[subroutine_name]["inlined_lines"] = {}
        self.func_info[subroutine_name]["inlined_lines"][filename] = new_lines

def get_used_files(make_output,directory):
    files = []
    if make_output is not None:
        with open(make_output, mode="r") as file:
                lines = file.readlines()
                for line in lines:
                    search = re.search("([^\s]+\.f90)", line)
                    if search is not None:
                        files.append(f"{directory}/{search.group(1)}")
        return files
    return glob.glob(f"{directory}/**/*.f90",recursive=True)

def main():
    argparser = argparse.ArgumentParser(description="Tool to find static writes in Fortran code",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    argparser.add_argument("-f", "--function", help="function to be parsed", required=True)
    argparser.add_argument("-F", "--file", help="File where to parse the function", required=True)
    argparser.add_argument("-m", "--use-make-output",type=str, help="Pass a file path if want to only analyze the files used in build process. Takes in the print output of make. If not given traverses the file structure to find all .f90 files in /src")
    argparser.add_argument("-c", "--communication",default=True,action=argparse.BooleanOptionalAction, help="Whether to check for mpi calls, can be used to check if code to be multithreaded is safe")
    argparser.add_argument("-o", "--offload",default=False,action=argparse.BooleanOptionalAction, help="Whether to offload the current code. Cannot multithread and offload at the same time")
    argparser.add_argument("-t", "--test",default=False,action=argparse.BooleanOptionalAction, help="Whether to generate an inlined version of the given subroutine (mainly used for testing)")
    argparser.add_argument("-d", "--directory",required=True, help="From which directory to look for files")
    argparser.add_argument("--sample-dir", required=False, help="Sample")
    # argparser.add_argument("-od", "--out-directory",required=True, help="To which directory write include files")
    argparser.add_argument("-M", "--Makefile", help="Makefile.local from which used modules are parsed")
    argparser.add_argument("-b", "--boundcond",default=False,action=argparse.BooleanOptionalAction, help="Whether the subroutine to offload is a boundcond or not")
    argparser.add_argument("-s", "--stencil",default=False,action=argparse.BooleanOptionalAction, help="Whether the subroutine to offload is a stencil op e.g. RK3 or not")
    argparser.add_argument("--to-c",default=False,action=argparse.BooleanOptionalAction, help="Whether to translate offloadable function to single threaded C for testing")
    
    args = argparser.parse_args()
    config = vars(args)
    filename = config["file"]
    directory_name = config["directory"]
    files = get_used_files(config["use_make_output"],directory_name)
    files = [file for file in files if os.path.isfile(file)]
    parser = Parser(files,config)
    # lines = [line.replace("\n","").strip() for line in open("res.txt").readlines()]
    # lines = parser.unroll_constant_loops(lines,{})
    # print("\n\nres\n\n")
    # for line in lines:
    #     print(line)
    # exit()
    # lines = parser.eliminate_while(lines)
    # for line in lines:
    #     print(line)
    # exit()
    # lines = parser.get_lines(file)
    # for line in lines:
    #     print(line)
    # print(parser.func_info["calc_pencils_magn_mf"]["files"])
    # exit()
    # line = "(.true. .or. .false.) .and. nwgrid/=1"
    # print(line)
    # print("-->")
    # print(self.evaluate_boolean(line))
    # exit()

    # lines = [line.replace("\n","") for line in open("./res.txt","r").readlines()]
    # lines = parser.transform_case(lines)
    # lines = parser.replace_var_in_lines(lines,"mass_source_profile","'exponential'")
    # local_variables = {}
    # parser.evaluate_ifs(lines,local_variables)
    # lines= parser.eliminate_while(lines)
    # print("res lines")
    # for line in lines:
    #     print(line)
    # exit()

    for module in parser.module_info:
        parser.parse_module(module)

    variables = {}
    # parser.get_lines(filename)
    # parser.get_lines("./test-dir/pencil-code/src/noyinyang.f90")
    # print(parser.func_info["in_overlap_mask"]["lines"].keys())
    # exit()
    header = "./pencil-code/src/mpicomm.h"
    subroutine_name = config["function"]
    # cparam_pencils.inc might not be included
    parser.load_static_variables(f"{parser.directory}/cparam.f90")
    #PENCIL SPECIFIC
    parser.populate_pencil_case()
    parser.static_variables["lpencil"] = {"type": "logical", "dims": ["npencils"], "threadprivate": False}
    parser.static_variables["nghost"] = {"type": "integer", "dims": [],"threadprivate": False}
    # if "lpencil" not in parser.static_variables:
    #     parser.load_static_variables("./backup/pencil-code/samples/gputest/src/cparam.inc")
    #     parser.load_static_variables("./pencil-code/src/cparam_pencils.inc")
    #     parser.load_static_variables("./pencil-code/src/special/axionSU2back.f90")
    # if "lpencil" not in parser.static_variables:
    #     exit()

    # print("density_scale" in parser.static_variables)
    # if "density_scale" not in parser.static_variables:
    #     print("density_scale" not in parser.static_variables)
    #     exit()
    # else:
    #     print(parser.static_variables["density_scale"])
    if config["test"] or config["offload"]:
        lines = parser.get_lines(f"{parser.sample_dir}/src/cparam.inc")
        lines = parser.get_lines(f"{parser.sample_dir}/src/cparam.inc")
        writes = parser.get_writes(lines,False)
        for write in writes:
            if write["variable"][0] == "l" and (write["value"] == ".true." or write["value"] == ".false."):
                parser.flag_mappings[write["variable"]] = write["value"]
            elif write["variable"] == "nghost":
                parser.flag_mappings[write["variable"]] = write["value"]

        
        #get flag mappings from cdata
        lines = parser.get_lines(f"{parser.sample_dir}/src/cdata.f90")
        writes = parser.get_writes(lines,False)
        for write in writes:
            if write["variable"][0] == "l":
                if write["variable"] not in parser.default_mappings:
                    parser.default_mappings[write["variable"]] = write["value"]

        #get flags from params.log
        lines = parser.get_lines(f"{parser.sample_dir}/data/params.log")

        for module in ["grav_","density_","magnetic_","hydro_"]:
            for prefix in ["run"]:
                grav_lines = parser.get_pars_lines(prefix, module,lines)
                res_lines = []
                for line in grav_lines:
                    res_lines.extend([part.strip() for part in line[0].split(",")])
                grav_lines = [(line,0) for line in res_lines]
                grav_writes = parser.get_writes(grav_lines,False)
                for write in grav_writes:
                    if write["value"] == "t":
                        parser.flag_mappings[write["variable"]] = ".true."
                    elif write["value"] == "f":
                        parser.flag_mappings[write["variable"]] = ".false."
                    elif "'" in write["value"]:
                        parsed_value = "'" + write["value"].replace("'","").strip() +"'"
                        parser.flag_mappings[write["variable"]] = parsed_value
                    #impossible value
                    elif write["value"] == "3.908499939e+37":
                        parser.flag_mappings[write["variable"]] = "impossible"
                    else:
                        parser.flag_mappings[write["variable"]] = write["value"]
        for map_param in parser.default_mappings:
            if map_param not in parser.flag_mappings:
                parser.flag_mappings[map_param] = parser.default_mappings[map_param]

        #for testing
        # parser.flag_mappings["topbot"] = "top"
        # parser.flag_mappings["lone_sided"] = ".false."
        # parser.flag_mappings["loptest(lone_sided)"] = ".false."

        if config["stencil"]:
            parser.offload_type = "stencil"
            #der and others are handled by the DSL
            parser.ignored_subroutines.extend(["der","der6"])
            #for testing excluded
            parser.ignored_subroutines.extend(["calc_slope_diff_flux"])
            parser.safe_subs_to_remove.extend(["calc_slope_diff_flux"])
        elif config["boundcond"]:
            parser.offload_type = "boundcond"

    if config["test"]:
        parser.inline_all_function_calls(filename,subroutine_name,add_init_lines=True) 
        new_lines = parser.func_info[subroutine_name]["inlined_lines"][filename]
        local_variables = {parameter:v for parameter,v in parser.get_variables(new_lines, {},filename).items() }
        writes = parser.get_writes(new_lines)
        #adding expand size in test
        new_lines = [parser.expand_size_in_line(line[0],local_variables,writes) for line in new_lines]
        new_lines= parser.eliminate_while(new_lines)
        new_lines = parser.unroll_constant_loops(new_lines,local_variables)
        global_loop_lines,iterators = parser.get_global_loop_lines(new_lines,local_variables)
        if len(global_loop_lines) > 0:
            new_lines= parser.remove_global_loops(new_lines,local_variables,global_loop_lines,iterators)
        # new_lines= parser.inline_1d_writes(new_lines,local_variables)
        print("\n\nDONE transform\n\n")
        print("LEN GLOBAL LOOP LINE",len(global_loop_lines))

        out_file = open(f"{parser.directory}/inlined_bc.inc","w")
        for line in new_lines:
            res_line = line.replace(subroutine_name,f"{subroutine_name}_inlined")
            print(res_line)
            out_file.write(f"{res_line}\n")
        out_file.close()

        out_file = open(f"{parser.directory}/inlined_bc_declaration.inc","w")
        out_file.write("public:: " +  subroutine_name + "_inlined")
        out_file.close()

        #generate new boundcond.f90
        old_lines = open(f"{parser.directory}/boundcond.f90").readlines()
        res_lines = []
        for line_index,line in enumerate(old_lines):
            res_line = line.strip()
            #don't want to take line continuation into account
            if len(res_line) > 0 and res_line[0] != "!" and subroutine_name in res_line and "&" not in res_line:
                func_calls = [call for call in parser.get_function_calls_in_line(res_line,{}) if call["function_name"] == subroutine_name]
                if len(func_calls) == 1:
                    print("should replace func call in line", line)
                    print("with test line")
                    print("--->")
                    replace_val = f"test_bc({subroutine_name},{subroutine_name}_inlined"
                    for param in func_calls[0]["parameters"]:
                        replace_val += "," + param
                    replace_val += ")"
                    res_line = parser.replace_func_call(res_line,func_calls[0],replace_val)
                    print(res_line)
            # if "endmodule" in res_line:
            #     res_lines.append("include 'inlined_bc'")
            res_lines.append(res_line)
        # out_file = open(f"{parser.directory}/boundcond.f90","w")
        # for line in res_lines:
        #     out_file.write(f"{line}\n")
        # out_file.close()
        print("done test setup")
        exit()
    if config["offload"]:
        #get flag mappings from cparam.inc
        #for testing
        # parser.flag_mappings["topbot"] = "top"
        # parser.flag_mappings["lone_sided"] = ".true."
        #
        parser.ignored_subroutines.extend(["mpiwtime","random_number_wrapper"])
        parser.ignored_subroutines.extend(["sum_mn_name","max_mn_name","yzsum_mn_name_x","xzsum_mn_name_y","xysum_mn_name_z","zsum_mn_name_xy","ysum_mn_name_xz","phizsum_mn_name_r"])
        parser.safe_subs_to_remove.extend(["sum_mn_name","max_mn_name","yzsum_mn_name_x","xzsum_mn_name_y","xysum_mn_name_z","zsum_mn_name_xy","ysum_mn_name_xz","phizsum_mn_name_r"])
        parser.ignored_subroutines.extend(["diagnostic_magnetic","xyaverages_magnetic","yzaverages_magnetic","xzaverages_magnetic"])
        parser.safe_subs_to_remove.extend(["diagnostic_magnetic","xyaverages_magnetic","yzaverages_magnetic","xzaverages_magnetic"])
        parser.ignored_subroutines.extend(["timing"])
        parser.safe_subs_to_remove.extend(["timing"])
        parser.ignored_subroutines.extend(["calc_diagnostics_density","calc_diagnostics_magnetic","calc_diagnostics_energy"])
        parser.safe_subs_to_remove.extend(["calc_diagnostics_density","calc_diagnostics_magnetic","calc_diagnostics_energy"])

        parser.ignored_subroutines.extend(["die_gracefully", "stop_it","stop_it_if_any"])
        parser.safe_subs_to_remove.extend(["die_gracefully","stop_it","stop_it_if_any"])
        parser.safe_subs_to_remove.extend(["open","close"])

        parser.ignored_subroutines.extend(["fatal_error","not_implemented","fatal_error_local","error","inevitably_fatal_error"])
        parser.safe_subs_to_remove.extend(["fatal_error","not_implemented","fatal_error_local","error","inevitably_fatal_error"])

        parser.safe_subs_to_remove.extend(["write"])

        parser.inline_all_function_calls(filename,subroutine_name,add_init_lines=True) 
        new_lines = parser.func_info[subroutine_name]["inlined_lines"][filename]
        print("\n\nDONE inlining\n\n")
        local_variables = {parameter:v for parameter,v in parser.get_variables(new_lines, {},filename).items() }
        transform_func = parser.transform_line_boundcond if config["boundcond"] else parser.transform_line_stencil
        res = parser.transform_lines([line[0] for line in new_lines],local_variables,transform_func)
        file = open("res.cuh","w")
        for line in res:
            file.write(f"{line}\n") 
        file.close()
        print("DONE")

        print("Deduce dims")
        print("Ranges:")
        has_x_range = False
        has_y_range = False
        has_z_range = False
        for range,index,line in parser.ranges:
            if range not in [":","1:mx","1:my","1:mz"]:
                print(range)
                print("NOT full range check what to do?")
                exit()
            has_x_range = has_x_range or (index == 0)
            has_y_range = has_y_range or (index == 1)
            has_z_range = has_z_range or (index == 2)
        x_dim= "mx" if has_x_range else "1"
        y_dim= "my" if has_y_range else "1"
        z_dim= "mz" if has_z_range else "1"
        print(f"{x_dim},{y_dim},{z_dim}")
        print(config)
        print(parser.test_to_c)
        exit()

    check_functions = []
    if config["communication"]:
        check_functions = parser.get_function_declarations(parser.get_lines(header), [], header)

    # parser.parse_subroutine_all_files(subroutine_name, "", check_functions, False, {})

    #slice buffers are safe, TODO add check to proof this
    parser.ignored_subroutines.extend(["store_slices","store_slices_scal","store_slices_vec"])
    parser.parse_subroutine_in_file(filename, subroutine_name, check_functions, config["offload"])
    print("DONE PARSING")
    print("modules")
    for module in parser.module_variables:
        print(module)
    parser.save_static_writes(parser.static_writes)

    writes = parser.static_writes
    variables = list(set([x["variable"] for x in writes if not x["local"]]))
    variables = list(set(variables))

    #df and f are accessed in threadsafe way
    critical_variables = ["df","f","num_of_diag_iter_done"]

    public_variables = [variable for variable in variables if variable in parser.static_variables and parser.static_variables[variable]["public"]]
    threadprivate_declarations = parser.get_threadprivate_declarations_and_generate_threadpriv_modules(parser.used_files,variables,critical_variables)
    threadprivate_var = []

    all_file = open(f"{parser.directory}/omp_includes/omp.inc","w")
    for file in threadprivate_declarations:
        include_filename = f"{parser.directory}/omp_includes/{file.rsplit('/',1)[-1].replace('.f90','')}_omp.inc"
        print("Creating ", include_filename)
        include_file = open(include_filename, "w")
        include_file.write(threadprivate_declarations[file])
        all_file.write(threadprivate_declarations[file])
        include_file.close()
    all_file.close()

    print("DONE")
    exit()

            

    

if __name__ == "__main__":
    main()
