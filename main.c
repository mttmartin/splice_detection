#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PTHREADS 1
#ifdef PTHREADS
#include <pthread.h>
#endif

#define BUFFER_MAX 1024


struct ParamContainer
{
    int seed_len;
    
    int part1_small_size_threshold;
    int part2_small_size_threshold;
    int mM1_s_threshold;
    
    float mM1_ratio_roof;
    float mM2_ratio_floor;
    
    int last_n_length;
    float last_n_mM_threshold;
};

enum {GOT_LAST_SEQS, FERROR, CONTINUE};

#ifdef PTHREADS

struct pthread_args {
	char *genome;
	struct ParamContainer paramcontainer;
};

pthread_mutex_t queue_mutex = PTHREAD_MUTEX_INITIALIZER;


pthread_mutex_t finished_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t finished_cond = PTHREAD_COND_INITIALIZER;

pthread_mutex_t start_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t start_cond = PTHREAD_COND_INITIALIZER;


int count = 0;
int term = 0;
int finished = 0;

int max_queue_size=0;
char **read_queue;
#endif




int array_sum(int array[], int length)
{
    int i;
    int ret = 0;
    for (i=0; i < length; ++i)
        ret = ret + array[i];
    return ret;
}

int strncpy_to_newline(char *string1, char *string2, int len)
{
    int i = 0;
    
    while (*string2 != '\0' && *string2 != '\n' && i < len)
    {
        *string1 = *string2;
        string1++;
        string2++;
        i++;
    }
    *string1 = '\0';
    return 0;
}

int get_genome(char *filename, char **genome)
{
	FILE *f = fopen(filename, "rb");
	long file_size;

	if (f == NULL)
	{
		fprintf(stderr, "Error opening genome");
		free(*genome);
		return -1;
	}

	fseek(f, SEEK_SET, SEEK_END);
	file_size = ftell (f);	
	*genome = malloc(sizeof(char)*file_size);

	fseek(f, SEEK_CUR, SEEK_SET);
	fread(*genome, 1, file_size, f);


	fclose(f);
	return 0;
}


int is_valid_read(char *read)
{	
	char c;
	int i = 0;
	while ( (read[i] != '\0') && (read[i] != '\n'))
	{
		c = read[i];
		if ((c == 'A') || (c == 'T') || (c == 'G') || (c == 'C'))
		{
			i++;
		}
		else
			return 0;
	}	
	return 1;
}


int get_next_sequences(char *file_name, int number, char buffers[][1024], fpos_t *last_pos)
{
	FILE *f = fopen(file_name, "r");
    
	if (f == NULL)
	{
		fprintf(stderr, "Error opening sequence file\n");
		return FERROR;
	}
    
	fsetpos(f, last_pos);

	int i;
	char line_buffer[1024] = {'\0'};
	for (i=0; i <= number; i++)
	{
		if (fgets(line_buffer,1024,f) != NULL)
		{
			switch (line_buffer[0])
			{
				case '>':
				case ';':
				{
					i--;
					break;
				}
				
				default:
					strncpy_to_newline(buffers[i], line_buffer, strlen(line_buffer));
			}
        }
		else
		{
			while (i <= number)
			{
				*buffers[i] = '\0';
				i++;
			}
			return GOT_LAST_SEQS;
		}
	}

	fgetpos(f, last_pos);
	fclose(f);
	return CONTINUE;
}

void get_genome_chunk(int location, unsigned long seq_len, char *genome, char *return_buffer)
{
    strncpy(return_buffer, genome+location, seq_len);
}


int location_in_genome(char *seq, char *genome)
{
    int max_mismatch = 1;
    int i;
    
    for (i=0; i < strlen(genome) - strlen(seq); i++)
    {
        int mismatches = 0;
        int j;
        for (j=0; j < strlen(seq); j++)
        {
            if (seq[j] != genome[i+j])
                mismatches++;
                
        }
        if (mismatches < max_mismatch)
            return i;
    }
    
    return -1;
}

int is_in_genome(char *seq, char *genome)
{
    if (location_in_genome(seq, genome) >= 0)
        return 1;
    else
        return 0;
}

int determine_splice_loc(char *seq, char *chunk, struct ParamContainer paramcontainer)
{

    int part1_len = 0;
    int part2_len = 0;
    int mM1 = 0;
    int mM2 = 0;
    int mM1_s = 0;
    int part = 1;
    
    int last_n_list[paramcontainer.last_n_length];
    int last_n_iter = 0;
    
    int i;
    for (i=0; i <= paramcontainer.last_n_length; ++i)
    {
        last_n_list[i] = 0;
    }
    
    int splice_location = -1;

    

    for (i=0; i < strlen(seq); i++)
    {
        if (part == 1)
        {
            part1_len++;
            if (seq[i] == chunk[i])
            {
                mM1_s = 0;
                last_n_list[last_n_iter] = 0;
            }
            else
            {
                mM1_s++;
                mM1++;
                
                last_n_list[last_n_iter] = 1;
                
                if (part1_len < paramcontainer.seed_len)
                    return -1;
            }

            if ((part1_len >= paramcontainer.part1_small_size_threshold) && ((float)mM1/part1_len > paramcontainer.mM1_ratio_roof))
                return -1;
            
            if ((part1_len >= paramcontainer.part1_small_size_threshold) && ((mM1_s > paramcontainer.mM1_s_threshold) || (array_sum(last_n_list, paramcontainer.last_n_length)/ (float)paramcontainer.last_n_length > paramcontainer.last_n_mM_threshold)))
            {
                i = i - mM1_s;
				part = 2;
				part1_len = part1_len - mM1_s;
				splice_location = i+1;
            }
            
            
            if (last_n_iter < paramcontainer.last_n_length)
                last_n_iter = 0;
            else
                last_n_iter++;
            
        }
        else if (part == 2)
        {
            part2_len++;
            if (seq[i] != chunk[i])
                mM2++;
            
            if ((part2_len > paramcontainer.part2_small_size_threshold) && ( (float) mM2/part2_len < paramcontainer.mM2_ratio_floor))
                return -1;
            
        }
    }

    return splice_location;
}

void check_for_splice(char *read, char *genome, struct ParamContainer paramcontainer)
{
	//fprintf(stderr, "check_for_splice-->Seq:%s\n", read);
    unsigned long read_len = strlen(read);
   

	if (!is_valid_read(read))
		return;
    
	int i;
    for (i=0; i < strlen(genome)-read_len; i++)
    {
	    char *genome_chunk = calloc(read_len+1, sizeof(char));
		get_genome_chunk(i, read_len, genome, genome_chunk);
        int splice_site = determine_splice_loc(read, genome_chunk, paramcontainer);
        
        if (splice_site > 0)
        {

            int part1_loc = 0;
            int part2_loc = 0;
		
			
			char *part1 = calloc(read_len, sizeof(char));
			char *part2 = calloc(read_len, sizeof(char));

			// cut read into parts
            strncpy(part1, read, splice_site);
            
            int j;
            for (j=0; !(is_in_genome(part1, genome)) && ((int)splice_site-j >= 0); j++)
            {
                    memset(part1, 0, read_len+1);
					strncpy(part1, read, splice_site-j);
            }
            
            if (strlen(part1) < paramcontainer.part1_small_size_threshold)
            {
                free(genome_chunk);
                continue;
            }

            splice_site = (int) strlen(part1);
			
			strncpy(part2, read+splice_site, strlen(read)-splice_site);
			
			part1_loc = location_in_genome(part1, genome);
            part2_loc = location_in_genome(part2, genome);
            
            if ((is_in_genome(part2, genome)) && (part1_loc < part2_loc))
            {
				int read_pol = '-'; // Dummy value until polarity is added
				int donor_loc = part1_loc + strlen(part1);
				int acceptor_loc = part2_loc;
				int intron_len = acceptor_loc - donor_loc;
                
				printf("%s,%i,%i,%i,%lu,%i,%c\n", read, donor_loc, acceptor_loc, splice_site, read_len, intron_len, read_pol);
				
				free(genome_chunk);
                free(part1);
                free(part2);
                return;
            }

        free(part1);
		free(part2); 
        }
       
	   	free (genome_chunk);

    }
    
}

#ifdef PTHREADS
void *pthread_worker (void *data)
{
	
	struct pthread_args *p_args = (struct pthread_args *) data;

	pthread_mutex_lock(&start_mutex);
	pthread_cond_wait(&start_cond, &start_mutex);
	pthread_mutex_unlock(&start_mutex);

	while (1)
	{
		if (term == 1)
		{
			return NULL;
		}

		pthread_mutex_lock(&queue_mutex);
		if (count > max_queue_size)
		{
			//printf("finished...\n");
			count = 0;
			pthread_mutex_lock(&finished_mutex);
			finished=1;
			pthread_cond_broadcast(&finished_cond);
			pthread_mutex_unlock(&finished_mutex);
			pthread_mutex_lock(&start_mutex);
			//printf("waiting for start\n");
			pthread_cond_wait(&start_cond, &start_mutex);
			pthread_mutex_unlock(&start_mutex);
			pthread_mutex_unlock(&queue_mutex);
			finished=0;
			//printf("got start\n");
			continue;
		}

		// tmp created to ensure read is null terminated
		char *read = read_queue[count];
		if ((read == NULL) || (*read == '\0'))
		{
			pthread_mutex_unlock(&queue_mutex);
			count++;
			return NULL;
		}
		
		char *tmp = calloc(strlen(read)+1, sizeof(char));
		memcpy(tmp, read, strlen(read));
		count++;
		pthread_mutex_unlock (&queue_mutex);
		//printf("count:%i\n", count);

		check_for_splice (tmp, p_args->genome, p_args->paramcontainer);
		free(tmp);
	}

}
#endif



int main(int argc, char *argv[])
{
	double current_version = 0.1;
	char *file_loc = NULL;
	char *genome = NULL;
	char *status_log_file_loc = NULL;
	FILE *status_log_file = NULL;
	int log_interval = 1000;
	
	int thread_num = 5;
    struct ParamContainer paramcontainer;
    
    paramcontainer.seed_len = 3;
    paramcontainer.part1_small_size_threshold = 3;
    paramcontainer.part2_small_size_threshold = 3;
    paramcontainer.mM1_ratio_roof = 0.5;
    paramcontainer.mM2_ratio_floor = 0.2;
    paramcontainer.mM1_s_threshold = 3;
    
    paramcontainer.last_n_length = 4;
    paramcontainer.last_n_mM_threshold = 0.20;


	if (argc < 2)
	{
		fprintf(stderr, "Incorrect usage.\n");
		return 1;
	}

	int i;
	for (i=1; i < argc; i++)
	{
		// input file: -i, --input
		// //output file: -o, --output
		// genome file: -g, --genome 
		// print status log: -l, --status_log 
		// thread number: -t --threads
		// seed_len: -s --seed_len
		// part1_small_size_threshold: --partA_min_len
		// part2_small_size_threshold: --partB_min_len
		// mM1_ratio_roof: --mM1_roof
		// mM2_ratio_floor: --mM2_floor
		// mM1_s_threshold: --mMs
		// last_n_length: --last_n_len
		// last_n_mM_threshold: --mM_last_n

		if ( (strncmp(argv[i], "-", 1) == 0) || (strncmp(argv[i], "--", 2) == 0))
		{
			if ( (strcmp(argv[i], "--input") == 0) || (strcmp(argv[i], "-i") == 0))
			{
				file_loc = argv[i+1];
			}
			else if ( (strcmp(argv[i], "--genome") == 0) || (strcmp(argv[i], "-g") == 0))
			{
				get_genome(argv[i+1], &genome);
			}
			else if ( (strcmp(argv[i], "--status_log") == 0) || (strcmp(argv[i], "-l") == 0))
			{
				status_log_file_loc = malloc(sizeof(char)*strlen(argv[i+1]));
				strncpy(status_log_file_loc, argv[i+1], strlen(argv[i+1]));
			}
			else if ( (strcmp(argv[i], "--threads") == 0) || (strcmp(argv[i], "-t") == 0))
			{
				thread_num = atoi(argv[i+1]);
				fprintf(stderr, "threads:%s", argv[i+1]);
			}
			else if ( (strcmp(argv[i], "--seed_len") == 0) || (strcmp(argv[i], "-s") == 0))
			{
				paramcontainer.seed_len = atoi(argv[i+1]);
			}
			else if (strcmp(argv[i], "--mM1_roof") == 0)
			{
				if ((atof(argv[i+1]) > 0) && (atof(argv[i+1]) < 1))
					paramcontainer.mM1_ratio_roof = atof(argv[i+1]);
				else
					fprintf(stderr, "Warning: --mM1_roof must be between zero and one\n");
			}
			else if (strcmp(argv[i], "--mM2_floor") == 0)
			{
				if ((atof(argv[i+1]) > 0) && (atof(argv[i+1]) < 1))
					paramcontainer.mM2_ratio_floor = atof(argv[i+1]);
				else
					fprintf(stderr, "Warning: --mM2_floor must be between zero and one\n");
			}
			else if (strcmp(argv[i], "--mMs") == 0)
			{
 				paramcontainer.mM1_s_threshold = atoi(argv[i+1]);
			}
			else if (strcmp(argv[i], "--last_n_len") == 0)
			{
			    paramcontainer.last_n_length = atoi(argv[i+1]);
			}
			else if (strcmp(argv[i], "--mM_last_n") == 0)
			{
				if ((atof(argv[i+1]) > 0) && (atof(argv[i+1]) < 1))
					paramcontainer.last_n_mM_threshold = atof(argv[i+1]);
				else
					fprintf(stderr, "Warning: --mM_last_n must be between zero and one\n");
			}
			else if (strcmp(argv[i], "--partA_min_len") == 0)
			{
				paramcontainer.part1_small_size_threshold = 3;
			}
			else if (strcmp(argv[i], "--partB_min_len") == 0)
			{
				paramcontainer.part2_small_size_threshold = 3;
			}
			else if ((strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "--help") == 0))
			{
				fprintf(stderr, "Usage: splice_detection -i read_file -g genome_file [options]\n");
				fprintf(stderr, "See man splice_detection for further information\n");
				return 0;
			}
			else if ((strcmp(argv[i], "-v") == 0) || (strcmp(argv[i], "--version") == 0))
			{
				fprintf(stderr, "splice_detection %g\n", current_version);
				return 0;
			}
			else
			{
				fprintf(stderr, "Warning ignoring unknown command:%s\n", argv[i]);
			}
		}

	}
	
	if ((file_loc == NULL) || (genome == NULL))
	{
		fprintf(stderr, "Error: must specify input file and genome\n");
		return 1;
	}
	


#ifdef PTHREADS
	printf("read, donor_loc, acceptor_loc, splice_loc, read_len, intron_len, read_pol\n");

	struct timespec ts1, ts2;
	ts1.tv_sec = 0;
	ts1.tv_nsec = 500000L;
	max_queue_size = thread_num;
	struct pthread_args *args = malloc(sizeof(struct pthread_args));
	args->genome = genome;
	args->paramcontainer = paramcontainer;

	pthread_t thread_id[thread_num];
	read_queue = (char**)malloc(1+max_queue_size * sizeof(char*));
	
	for (i=0; i <= max_queue_size; ++i)
		read_queue[i] = malloc(sizeof(char)*1024);


	for (i=0; i < thread_num; i++)
	{
		//fprintf(stderr, "Creating thread:%i\n", i);
		pthread_create (&thread_id[i], NULL, pthread_worker,  (void*)args);
	}
	
	int j;
	int log_iter=0;
	int total_processed_reads=0;
	int get_next_sequences_ret;
	fpos_t last_pos = 0;
	char out_queue[max_queue_size][1024];
	
	int loop_num = 0;


	fprintf(stderr,"max_queue_size:%i\n", max_queue_size);
	while (1)
	{
		//fprintf(stderr,"starting loop:%i\n",loop_num);
		loop_num++;

		if (log_iter > log_interval)
		{
			if (status_log_file_loc)
			{
				status_log_file = fopen(status_log_file_loc, "w");
				fprintf(status_log_file, "%i\n", (loop_num*max_queue_size));
				fclose(status_log_file);
			}
			fprintf(stderr, "Log:%i\n", (loop_num*max_queue_size));
			log_iter = 0;
		}
		log_iter++;

		get_next_sequences_ret = get_next_sequences(file_loc, max_queue_size, out_queue, &last_pos);
		if (get_next_sequences_ret == FERROR)
		{
			fprintf(stderr, "Bailing due to file error\n");
			term=1;
			pthread_cond_broadcast(&start_cond);
			break;
		}

		for (j=0; j <= max_queue_size; ++j)
		{
			memset(read_queue[j], '\0', 1024);
			strncpy (read_queue[j], out_queue[j], strlen(out_queue[j]));
			//printf("read:%s\n", out_queue[j]);
		}
	
		pthread_cond_broadcast(&start_cond);
		nanosleep(&ts1, &ts2);
		pthread_mutex_lock(&finished_mutex);
		while (finished == 0)
			pthread_cond_wait(&finished_cond, &finished_mutex);

		pthread_mutex_unlock(&finished_mutex);
		if (get_next_sequences_ret == GOT_LAST_SEQS)
		{
			fprintf(stderr, "Got last batch of sequences\n");
			break;
		}
	}

	for (i=0; i < thread_num; i++)
	{
		//fprintf(stderr, "Joining thread:%i\n", i);
		pthread_join(thread_id[i], NULL);
	}
#endif

	
	free (read_queue);
	free(genome);

}
