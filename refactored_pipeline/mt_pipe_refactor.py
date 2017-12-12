#!/usr/bin/env python
# this pipeline is a mess.
# This program simply makes queue submissions to the cluster, and in doing so, eats a lot of down
import sys
import os
import os.path
import subprocess as sp
import multiprocessing
import mt_pipe_commands as mpcom
import mt_pipe_paths as mpfp
import time
Threads = str(multiprocessing.cpu_count())

run_jobs = False
def make_folder(folder_path):
    if not(os.path.exists(folder_path)):
        os.makedirs(folder_path)

class qsub_sync:
    #handles job syncing
    def check_job_finished(self, job_id, op_mode):
        if(op_mode == "completed"):
            jobs_list = sp.check_output(["showq", "-c"]).decode('ascii')
        elif(op_mode == "queue"):
            jobs_list = sp.check_output(["showq", "-i"]).decode('ascii')
        if (isinstance(job_id, int)):
            #single
            if(str(job_id) in jobs_list):
                return True
            else:
                return False
            
        elif(isinstance(job_id, list)):
            # multi
            match_count = 0
            for item in job_id:
                if(str(item) in jobs_list):
                    match_count += 1
            if(match_count != len(job_id)):
                return False
            else:
                return True
                
        elif(dependency_list is None):
            print("bad arg to check job finished")
            return False
        else:
            print("This isn't supposed to happen")
            sys.exit()
            
    def wait_for_sync(self, timeout, job_id, label):
        #does the actual wait for the qsub job to finish
        b_lock = True
        lockout_count = int(timeout)
        while(b_lock):
            if(self.check_job_finished(job_id, "completed")):
                b_lock = False
            else:
                time.sleep(1)
                #only start the countdown if the job is running, not while it's waiting in the queue.  
                #job could wait a long time in the queue if someone's clobbering the cluster
                if not(self.check_job_finished(job_id, "queue")):
                    lockout_count -= 1
            
            if(lockout_count <= 0):
                print(label, "took too long.  shutting down pipeline")
                sys.exit()

def main(input_folder, output_folder):
    # constants
    # -----------------------------
    single_mode = 0
    double_mode = 1

    # system vars
    # -------------------------------
    # operating mode:
    # 0: single-ended
    # 1: paired
    # 2: error
    operating_mode = 0
    sync_obj = qsub_sync()
    start_time = time.time()
    #note: this also needs to support paired and single-ended data
    #input folder is the main location of the dump.
    #
    #for genome in sorted(os.listdir(input_folder)):
    file_list = []
    Network_list = []
    #only seems to look for *1.fastq, and nothing else.  the whole loop is wasting time.  
    raw_sequence_path = input_folder + "/raw_sequences/"
    if not os.path.exists(raw_sequence_path):
        print("No sequences found.  to use pipeline, please place fastq file at:", raw_sequence_path)
        os.makedirs(raw_sequence_path)
        sys.exit()
    else:
        # folder found.  now see if it's single-ended or paired
        # for single-ended, have only 1 file.  if doubled, have both files.  Else, stop
        genome_file_count = len(os.listdir(raw_sequence_path))
        print("number of files:", genome_file_count)
        if(genome_file_count == 1):
            print("OPERATING IN SINGLE-ENDED MODE")
        elif(genome_file_count == 2):
            print("OPERATING IN PAIRED-MODE")
        else:
            print("Too many genome files here.  Get rid of all non-essentials")
            sys.exit()
        
        operating_mode = genome_file_count - 1
        
        # is the file too big?
        # split it.
        
        #init a command object, and start making commands
        #sys.exit()
        
        if(operating_mode == double_mode):
            
            preprocess_label = "preprocess"
            raw_pair_0_path = raw_sequence_path + sorted(os.listdir(raw_sequence_path))[0]
            raw_pair_1_path = raw_sequence_path + sorted(os.listdir(raw_sequence_path))[1]
            comm = mpcom.mt_pipe_commands(Quality_score = 33, Thread_count = 16, raw_sequence_path_0 = raw_pair_0_path, raw_sequence_path_1 = raw_pair_1_path)
            preprocess_job_id = comm.create_pbs_and_launch("preprocess", comm.create_pre_double_command(preprocess_label), run_job = True)
            
            rRNA_filter_job_id = []
            rRNA_filter_job_id.append(comm.create_pbs_and_launch("rRNA_filter", comm.create_rRNA_filter_prep_command("rRNA_filter", 5, "preprocess"), dependency_list = preprocess_job_id, run_job = True))
            
            #standalone
            #rRNA_filter_job_id.append(comm.create_pbs_and_launch("rRNA_filter", comm.create_rRNA_filter_prep_command("rRNA_filter", 5, "preprocess"), run_job = True))
            #--------------------------
            #print("-----------------------------------")
            #print(rRNA_filter_job_id)
            #print(rRNA_filter_job_id[0])
            #print("-----------------------------------")
            
            #constructing folder paths, so we have an easier time
            rRNA_filter_orphans_fasta_folder = os.getcwd() + "/rRNA_filter/data/orphans/orphans_fasta/"
            rRNA_filter_pair_1_fasta_folder = os.getcwd() + "/rRNA_filter/data/pair_1/pair_1_fasta/"
            rRNA_filter_pair_2_fasta_folder = os.getcwd() + "/rRNA_filter/data/pair_2/pair_2_fasta/"
            #-------------------------------------------------------------------
            rRNA_filter_orphans_fastq_folder = os.getcwd() + "/rRNA_filter/data/orphans/orphans_fastq/"
            rRNA_filter_pair_1_fastq_folder = os.getcwd() + "/rRNA_filter/data/pair_1/pair_1_fastq/"
            rRNA_filter_pair_2_fastq_folder = os.getcwd() + "/rRNA_filter/data/pair_2/pair_2_fastq/"
            #--------------------------------------------------------------------
            rRNA_filter_orphans_mRNA_folder = os.getcwd() + "/rRNA_filter/data/orphans/orphans_mRNA/"
            rRNA_filter_pair_1_mRNA_folder =  os.getcwd() + "/rRNA_filter/data/pair_1/pair_1_mRNA/"
            rRNA_filter_pair_2_mRNA_folder =  os.getcwd() + "/rRNA_filter/data/pair_2/pair_2_mRNA/"
            #---------------------------------------------------------------------
            rRNA_filter_orphans_rRNA_folder = os.getcwd() + "/rRNA_filter/data/orphans/orphans_rRNA/"
            rRNA_filter_pair_1_rRNA_folder =  os.getcwd() + "/rRNA_filter/data/pair_1/pair_1_rRNA/"
            rRNA_filter_pair_2_rRNA_folder =  os.getcwd() + "/rRNA_filter/data/pair_2/pair_2_rRNA/"
            
            rRNA_filter_final_mRNA_folder = os.getcwd() + "/rRNA_filter/data/final_result/mRNA/"
            rRNA_filter_final_rRNA_folder = os.getcwd() + "/rRNA_filter/data/final_result/rRNA/"
            
            make_folder(rRNA_filter_final_mRNA_folder)
            make_folder(rRNA_filter_final_rRNA_folder)
            
            #this is gonna be hacky.... we have to wait until the prep stage is finished, but there's no nice way to sense it, through qsub
            #why delay?  because the following code needs the files present to generate the correct job.  
            #this delay, and timeout are a check against super large runaway jobs.
            
            sync_obj.wait_for_sync(600, rRNA_filter_job_id[0], "rRNA filter")
            print("moving onto INFERNAL")
            
            for item in os.listdir(rRNA_filter_orphans_fastq_folder):
                file_root_name = item.split('.')[0]
                rRNA_filter_job_id.append(
                    comm.create_pbs_and_launch(
                        "rRNA_filter", 
                        comm.create_rRNA_filter_command("rRNA_filter", "orphans", file_root_name), 
                        inner_name = file_root_name + "_infernal",
                        #dependency_list = rRNA_filter_job_id[0],
                        run_job = True
                    )
                )
                
            for item in os.listdir(rRNA_filter_pair_1_fastq_folder):
                file_root_name = item.split('.')[0]
                rRNA_filter_job_id.append(
                    comm.create_pbs_and_launch(
                        "rRNA_filter", 
                        comm.create_rRNA_filter_command("rRNA_filter", "pair_1", file_root_name), 
                        inner_name = file_root_name + "_infernal",
                        #dependency_list = rRNA_filter_job_id[0],
                        run_job = True
                    )
                )
                
            for item in os.listdir(rRNA_filter_pair_2_fastq_folder):
                file_root_name = item.split('.')[0]
                rRNA_filter_job_id.append(
                    comm.create_pbs_and_launch(
                        "rRNA_filter", 
                        comm.create_rRNA_filter_command("rRNA_filter", "pair_2", file_root_name), 
                        inner_name = file_root_name + "_infernal",
                        #dependency_list = rRNA_filter_job_id[0],
                        run_job = True
                    )
                )
            #wait for infernal to finish running
            sync_obj.wait_for_sync(600, rRNA_filter_job_id, "rRNA filter")
            #then we need to combine the splits into per-category files
            #this shouldn't be a qsub job
            
            cat_orphans_mRNA = "cat " + rRNA_filter_orphans_mRNA_folder + "* 1>>" + rRNA_filter_final_mRNA_folder + "orphans.fastq"
            cat_orphans_rRNA = "cat " + rRNA_filter_orphans_rRNA_folder + "* 1>>" + rRNA_filter_final_rRNA_folder + "orphans.fastq"
            
            cat_pair_1_mRNA = "cat " + rRNA_filter_pair_1_mRNA_folder + "* 1>>" + rRNA_filter_final_mRNA_folder + "pair_1.fastq"
            cat_pair_1_rRNA = "cat " + rRNA_filter_pair_1_rRNA_folder + "* 1>>" + rRNA_filter_final_rRNA_folder + "pair_1.fastq"
            
            cat_pair_2_mRNA = "cat " + rRNA_filter_pair_2_mRNA_folder + "* 1>>" + rRNA_filter_final_mRNA_folder + "pair_2.fastq"
            cat_pair_2_rRNA = "cat " + rRNA_filter_pair_2_rRNA_folder + "* 1>>" + rRNA_filter_final_rRNA_folder + "pair_2.fastq"
            
            #--------------------------------------------------------------------------------------------------------
            # stage 3: 
            
            
            end_time = time.time()
            print("Total runtime:", end_time - start_time)
        elif(operating_mode == single_mode):
            print("not ready")

        """
            # Preprocessing
            

            create_pbs_job("Preprocess", Input_FName, COMMANDS_Pre)        
            if(run_jobs):
                JobID_Pre = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Preprocess.pbs"])

            
            create_pbs_job("rRNA_Submit", Input_FName, COMMANDS_rRNA)
            if(run_jobs):
                JobID_rRNA = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_rRNA_Submit.pbs", "-W", "depend=afterok:" + JobID_Pre.strip("\n")])


            create_pbs_job("Combine", Input_FName, COMMANDS_Combine)
            if(run_jobs):
                JobID_Combine = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Combine.pbs", "-W", "depend=afterok:" + JobID_rRNA.strip("\n")])

                subprocess.call(["qalter", "-v", "JOB2=" + JobID_Combine.strip("\n").split(".")[0], JobID_rRNA.strip("\n")])

            #BIGG Database, AGORA Nature paper, Additional functionality
            # Transcript Assembly
            Contigs = os.path.join(Input_Path, os.path.splitext(Input_FName)[0] + "_SpadesOut", "contigs.fasta")

            create_pbs_job("Assemble", Input_FName, COMMANDS_Assemble)        
            if(run_jobs):
                JobID_Assemble = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Assemble.pbs", "-W", "depend=afterok:" + JobID_Combine.strip("\n")])

            # Protein Annotation BWA
            create_pbs_job("Annotate_BWA", Input_FName, COMMANDS_Annotate_BWA, "med")        
            if(run_jobs):
                JobID_Annotate_BWA = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_BWA.pbs", "-W", "depend=afterok:" + JobID_Assemble.strip("\n")])

            # Protein Annotation BLAT 1
            create_pbs_job("Annotate_BLAT1", Input_FName, COMMANDS_Annotate_BLAT1, "med")        
            if(run_jobs):
                JobID_Annotate_BLAT1 = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_BLAT1.pbs", "-W", "depend=afterok:" + JobID_Annotate_BWA.strip("\n")])

            # Protein Annotation BLAT 2
            create_pbs_job("Annotate_BLAT2", Input_FName, COMMANDS_Annotate_BLAT2, "med")
            if(run_jobs):
                JobID_Annotate_BLAT2 = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_BLAT2.pbs", "-W", "depend=afterok:" + JobID_Annotate_BWA.strip("\n")])

            # Protein Annotation BLAT 3
            create_pbs_job("Annotate_BLAT3", Input_FName, COMMANDS_Annotate_BLAT3, "med")       
            if(run_jobs):
                JobID_Annotate_BLAT3 = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_BLAT3.pbs", "-W", "depend=afterok:" + JobID_Annotate_BWA.strip("\n")])

            # Protein Annotation BLAT 4
            create_pbs_job("Annotate_BLAT4", Input_FName, COMMANDS_Annotate_BLAT4, "med")
            if(run_jobs):
                JobID_Annotate_BLAT4 = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_BLAT4.pbs", "-W", "depend=afterok:" + JobID_Annotate_BWA.strip("\n")])

            # Protein Annotation BLAT Postprocessing
            create_pbs_job("Annotate_BLAT_Postprocessing", Input_FName, COMMANDS_Annotate_BLAT_Post)
            if(run_jobs):
                JobID_Annotate_BLAT_Post = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_BLAT_Postprocessing.pbs", "-W", "depend=afterok:" + JobID_Annotate_BLAT1.strip("\n") + ":" + JobID_Annotate_BLAT2.strip("\n") + ":" + JobID_Annotate_BLAT3.strip("\n") + ":" + JobID_Annotate_BLAT4.strip("\n")])


            # Protein Annotation Diamond 1
            create_pbs_job("Annotate_DMD", Input_FName, COMMANDS_Annotate_Diamond1)        
            if(run_jobs):
                JobID_Annotate_Diamond1 = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_DMD1.pbs", "-W", "depend=afterok:" + JobID_Annotate_BLAT_Post.strip("\n")])

            # Protein Annotation Diamond 2
            create_pbs_job("Annotate_DMD2", Input_FName, COMMANDS_Annotate_Diamond2)  
            if(run_jobs):
                JobID_Annotate_Diamond2 = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_DMD2.pbs", "-W", "depend=afterok:" + JobID_Annotate_BLAT_Post.strip("\n")])

            # Protein Annotation Diamond 3
            create_pbs_job("Annotate_DMD3", Input_FName, COMMANDS_Annotate_Diamond3)
            if(run_jobs):
                JobID_Annotate_Diamond3 = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_DMD3.pbs", "-W", "depend=afterok:" + JobID_Annotate_BLAT_Post.strip("\n")])

            # Protein Annotation Diamond 4
          
            create_pbs_job("Annotate_DMD4", Input_FName, COMMANDS_Annotate_Diamond4)                
            if(run_jobs):
                JobID_Annotate_Diamond4 = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_DMD4.pbs", "-W", "depend=afterok:" + JobID_Annotate_BLAT_Post.strip("\n")])

            # Protein Annotation Diamond Postprocess
            create_pbs_job("Annotate_DMD_Postprocessing", Input_FName, COMMANDS_Annotate_Diamond_Post)
            if(run_jobs):
                JobID_Annotate_Diamond_Post = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_DMD_Postprocess.pbs", "-W", "depend=afterany:" + JobID_Annotate_Diamond1.strip("\n") + ":" + JobID_Annotate_Diamond2.strip("\n") + ":" + JobID_Annotate_Diamond3.strip("\n") + ":" + JobID_Annotate_Diamond4.strip("\n")])
            
            # Classify Reads
            create_pbs_job("Classify", Input_FName, COMMANDS_Classify, "high")
            if(run_jobs):
                JobID_Classify = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Classify.pbs", "-W", "depend=afterok:" + JobID_Annotate_Diamond_Post.strip("\n")])
            
            # Prepare EC annotation files
            EC_Split = os.path.join(Input_Filepath + "_EC_Annotation", "Split")
            EC_Output = os.path.join(Input_Filepath + "_EC_Annotation", "Output")
            
            create_pbs_job("EC_Preprocess", Input_FName, COMMANDS_EC_Preprocess)
            if(run_jobs):
                JobID_EC_Preprocess = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_EC_Preprocess.pbs", "-W", "depend=afterok:" + JobID_Annotate_Diamond_Post.strip("\n")])

            #EC detection
            create_pbs_job("Detect", Input_FName, COMMANDS_Detect)        
            if(run_jobs):
                JobID_Detect = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Detect.pbs", "-W", "depend=afterok:" + JobID_EC_Preprocess.strip("\n")])

            #Combine the EC files
            create_pbs_job("Combine_Detect", Input_FName, COMMANDS_Combine_Detect)
            if(run_jobs):
                JobID_Combine_Detect = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Combine_Detect.pbs", "-W", "depend=afterok:" + JobID_Detect.strip("\n")])
            
                subprocess.call(["qalter", "-v", "JOB2=" + JobID_Combine_Detect.strip("\n").split(".")[0], JobID_Detect.strip("\n")])

            
            #PRIAM stage
            create_pbs_job("PRIAM", Input_FName, COMMANDS_PRIAM)        
            if(run_jobs):
                JobID_PRIAM = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_PRIAM.pbs", "-W", "depend=afterok:" + JobID_EC_Preprocess.strip("\n")])

            # EC Diamond
            create_pbs_job("EC_Diamond", Input_FName, COMMANDS_EC_Diamond)        
            if(run_jobs):
                JobID_EC_Diamond = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_EC_Diamond.pbs", "-W", "depend=afterok:" + JobID_EC_Preprocess.strip("\n")])

            # EC Annotation Compile
            

            with open(os.path.splitext(Input_FName)[0] + "_EC_Postprocess.pbs", "w") as PBS_script_out:
                for line in PBS_Submit_LowMem.splitlines():
                    if "NAME" in line:
                        line = line.replace("NAME", os.path.splitext(Input_FName)[0] + "_EC_Postprocess")
                    if "ERROR" in line:
                        line = line.replace("ERROR", os.path.splitext(Input_FName)[0] + "_EC_Postprocess_ERR")
                    if "OUTPUT" in line:
                        line = line.replace("OUTPUT", os.path.splitext(Input_FName)[0] + "_EC_Postprocess_OUT")
                    if "COMMANDS" in line:
                        PBS_script_out.write("\n".join(COMMANDS_EC_Postprocess))
                        break
                    PBS_script_out.write(line + "\n")
            if(run_jobs):
                JobID_EC_Postprocess = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_EC_Postprocess.pbs", "-W", "depend=afterok:" + JobID_Combine_Detect.strip("\n") + ":" + JobID_PRIAM.strip("\n") + ":" + JobID_EC_Diamond.strip("\n")])

            # Network Generation
            

            with open(os.path.splitext(Input_FName)[0] + "_Network.pbs", "w") as PBS_script_out:
                for line in PBS_Submit_LowMem.splitlines():
                    if "NAME" in line:
                        line = line.replace("NAME", os.path.splitext(Input_FName)[0] + "_Network")
                    if "ERROR" in line:
                        line = line.replace("ERROR", os.path.splitext(Input_FName)[0] + "_Network_ERR")
                    if "OUTPUT" in line:
                        line = line.replace("OUTPUT", os.path.splitext(Input_FName)[0] + "_Network_OUT")
                    if "COMMANDS" in line:
                        PBS_script_out.write("\n".join(COMMANDS_Network))
                        break
                    PBS_script_out.write(line + "\n")
            
            if(run_jobs):
                JobID_Network = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Network.pbs", "-W", "depend=afterok:" + JobID_EC_Postprocess.strip("\n") + ":" + JobID_Classify.strip("\n")])
            Network_list.append(JobID_Network.strip("\n"))

    if len(file_list) > 1:
        os.chdir(output_folder)
        Input_Filepath = os.path.join(output_folder, os.path.splitext(os.path.basename(file_list[0]))[0])[:-11]
        
        with open(os.path.splitext(os.path.basename(Input_Filepath))[0] + "_Join.pbs", "w") as PBS_script_out:
            for line in PBS_Submit_LowMem.splitlines():
                if "NAME" in line:
                    line = line.replace("NAME", os.path.splitext(os.path.basename(Input_Filepath))[0] + "_Join")
                if "ERROR" in line:
                    line = line.replace("ERROR", os.path.splitext(os.path.basename(Input_Filepath))[0] + "_Join_ERR")
                if "OUTPUT" in line:
                    line = line.replace("OUTPUT", os.path.splitext(os.path.basename(Input_Filepath))[0] + "_Join_OUT")
                if "COMMANDS" in line:
                    PBS_script_out.write("\n".join(COMMANDS_Join))
                    break
                PBS_script_out.write(line + "\n")
        
        if(run_jobs):
            JobID_Join = subprocess.check_output(["qsub", os.path.splitext(os.path.basename(Input_Filepath))[0] + "_Join.pbs", "-W", "depend=afterok:" + ":".join(Network_list)])
"""                

if __name__ == "__main__":
    input_folder = sys.argv[1]
    output_folder = sys.argv[2]
    #scinet_user_name = sys.argv[3]
    os.chdir(output_folder)
    main(input_folder, output_folder)