all_lines = []
for energy_index in [4,5,7,8,9]:
    with open(f"f{energy_index}.txt", "r") as file:
        count_line = 0
        add_line = ""
        info_lines = []

        collect_variables = set()
        collect_fit_types = set()

        #list_of_fit_types = ['Vv_Wd', 'Vv_Wd_scale', 'Vv_rv_Wd_rwd', 'Vv_rv_Wd_rwd_scale', 'Vv_rv_Wd_rwd_Vso_rso', 'Vv_rv_Wd_rwd_Vso_rso_scale', 'Vv_rv_av_Wd_rwd_awd', 'Vv_rv_av_Wd_rwd_awd_scale',]
        list_of_fit_types = ['Vv_Wd', 'Vv_rv_Wd_rwd', 'Vv_rv_av_Wd_rwd_awd', 'Vv_rv_Wd_rwd_Vso_rso', 'Vv_Wd_scale', 'Vv_rv_Wd_rwd_scale', 'Vv_rv_av_Wd_rwd_awd_scale', 'Vv_rv_Wd_rwd_Vso_rso_scale',]
        list_of_variables = ['scale', 'Vv', 'Wd', 'rv','rwd', 'awd', 'Vso', 'rso', 'av',]

        list_of_values = [None]*10
        list_of_errors = [None]*10

        for line in file:

            if line.startswith("E="):

                if len(add_line)>0:
                    for i, var in enumerate(list_of_variables):
                        if list_of_values[i]==None:
                            #add_line = add_line + f"[{var}] ,,   \n"
                            add_line = add_line + ",,   "
                        else:
                            #add_line = add_line + f"[{var}] " + list_of_values[i] + ", " + list_of_errors[i] + ",   \n"
                            add_line = add_line + list_of_values[i] + ", " + list_of_errors[i] + ",   "

                    info_lines.append(add_line)
                    all_lines.append(add_line)
                    add_line = ""
                    list_of_values = [None]*10
                    list_of_errors = [None]*10

                line_new = line[2:]
                stripped_line = line_new.lstrip()
                tokens = stripped_line.split()
                #add_line = add_line + tokens[0] + ", " + tokens[2] + ", "
                add_line = add_line + tokens[0] + ", " + f"{list_of_fit_types.index(tokens[2])}" + ", "
                collect_fit_types.add(tokens[2])

            elif line.startswith("# Var"):
                line_new = line[8:]
                stripped_line = line_new.lstrip()
                tokens = stripped_line.split()
                collect_variables.add(tokens[0])
                #add_line = add_line + tokens[0] + ", " + tokens[2][:-1] + ", " + tokens[6][:-1] + ", "
                v_index = list_of_variables.index(tokens[0])
                list_of_values[v_index] = tokens[2][:-1]
                list_of_errors[v_index] = tokens[6][:-1]

            elif line.startswith("# ChiSq/N"):
                line_new = line[11:]
                stripped_line = line_new.lstrip()
                tokens = stripped_line.split()
                add_line = add_line + tokens[0] + ", "

            elif line.startswith("# sigma_R"):
                line_new = line[11:]
                stripped_line = line_new.lstrip()
                tokens = stripped_line.split()
                add_line = add_line + tokens[0] + ", "
                #print(tokens)

            elif line.startswith("exp"):
                scale = line[4:].strip()
                add_line = add_line + scale + ", "

            count_line = count_line + 1
            #if count_line>50: break

        if len(add_line)>0:
            for i, var in enumerate(list_of_variables):
                if list_of_values[i]==None:
                    #add_line = add_line + f"[{var}] ,,   \n"
                    add_line = add_line + ",,   "
                else:
                    #add_line = add_line + f"[{var}] " + list_of_values[i] + ", " + list_of_errors[i] + ",   \n"
                    add_line = add_line + list_of_values[i] + ", " + list_of_errors[i] + ",   "

            info_lines.append(add_line)
            all_lines.append(add_line)
            add_line = ""
            list_of_values = [None]*10
            list_of_errors = [None]*10


    content = '\n'.join(info_lines)
    #print(content)
    with open(f"strip_f{energy_index}.csv", "w") as file:
        file.write(content)

all_content = '\n'.join(all_lines)
#print(all_content)
file_name = f"strip_all.csv"
with open(file_name, "w") as file:
    file.write(all_content)
print(file_name)

#print("list_of_fit_types =",list(collect_variables))
#print("list_of_variables =",list(collect_fit_types))
