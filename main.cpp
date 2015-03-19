#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <thread>
#include <mutex>
#include <sstream>
#include <string.h>
#include <seqan/find.h>
#define current_version "0.1.1"


using namespace std;


mutex cout_mutex;

enum {ON_DEFAULT, ON_REVERSE, ON_ANTI};

struct ORF_Cordinate
{
	int begin;
	int end;
	string id;
	
	ORF_Cordinate(int b, int e)
	{
		begin = b;
		end = e;
	}

	ORF_Cordinate(int b, int e, string s)
	{
		begin = b;
		end = e;
		id = s;
	}
};

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

	int intron_min;
	int intron_max;
	
	int do_anti_sense;
	int do_reverse;

	bool ignore_noncoding;
	bool circular_genome;

	vector<ORF_Cordinate> cords; 
};

class ThreadPool
{
	bool finished;
	vector<thread> threads;
	vector<function<void()>> work_queue;

	mutex queue_mutex;


	void worker();

	public:
		ThreadPool();
		ThreadPool(int threads_requested);
		void terminate();
		void add_work(function<void()> work);
	


};

ThreadPool::ThreadPool()
{
	int thread_num = thread::hardware_concurrency();
	for (int i=0; i < thread_num; ++i)
	{
		threads.push_back(thread(&ThreadPool::worker, this));
	}
	
}

ThreadPool::ThreadPool(int threads_requested)
{
	int thread_num;
	
	if (threads_requested > threads.max_size())
	{
		cout << "Error: Thread number is too large.\n";
		exit(1);
	}
	else if (threads_requested < 0)
	{
		cerr << "Error: Thread number must be greater than 0.\n";
		exit(1);
	}
	
	if (threads_requested == 0)
		thread_num = thread::hardware_concurrency();
	else
		thread_num = threads_requested;

	for (int i=0; i < thread_num; ++i)
	{
		threads.push_back(thread(&ThreadPool::worker, this));		
	}

}


void ThreadPool::terminate()
{
	finished = true;
	for (int i=0; i < threads.size(); ++i)
	{
		threads[i].join();
	}
}


void ThreadPool::add_work(function<void()> work)
{
	queue_mutex.lock();
	work_queue.push_back(work);
	queue_mutex.unlock();
}

	
void ThreadPool::worker()
{
	while((!finished) || (work_queue.size() > 0))
	{
		function<void()> work;
		
		queue_mutex.lock();
			
		if (work_queue.size() > 0)
		{
			work = work_queue.back();
			work_queue.pop_back();
		}

		queue_mutex.unlock();	

		if (work != NULL)
		{
			work();
		}
		// Our queue was empty... wait a bit
		else
		{
			this_thread::yield();	
		}
	}	
}

bool is_valid_read (string read)
{
	vector<char> allowed {'A', 'T', 'G', 'C'};
	
	for (int i=0; i < read.length(); i++)
	{
		if (find(allowed.begin(), allowed.end(), read[i]) != allowed.end())
			continue;
		// For now just ignore Ns; do not display error message
		else if (read[i] == 'N')
		{
			return false;
		}
		else
		{
			cerr << "Warning: Invalid read detected: " << read << endl;
			return false;
		}
	}

	return true;
}

/* get_genome(char *filename, char **genome)
 * Reads file given by filename and stores the nucleotide sequence in **genome.
 * Returns 0 if successful, otherwise returns < 0.
 */

int get_genome(string filename, string *genome)
{
	ifstream f;
	string line_buffer;
	
	f.open(filename);
	
	if (!f)
	{
		cerr << "Error opening genome file\n";
		exit(1);
	}

	while ( (getline(f, line_buffer)))
	{
		if (line_buffer.front() == '>')
			continue;
		
		if (!is_valid_read(line_buffer))
		{
			cerr << "Genome file contains illegal characters or is invalid.\n";
			exit(1);
		}
		
		genome->append(line_buffer);

	}

	f.close();
	return 0;
}


/* get_ORFs(string filename, Paramcontainer paramcontainer)
 * Retrieves ORF locations from GFF file at filename and stores
 * coordinates for all ORFs in paramcontainer for later use */
bool get_ORFS(string filename, ParamContainer *paramcontainer)
{
	string line_buffer;
	string field_buffer;
	stringstream stream;
	ifstream f;
	f.open(filename);
	
	if (!f)
	{
		cerr << "Error opening GFF file\n";
		exit(1);
	}

	while ( (getline(f, line_buffer)))
	{
		
		
		stream << line_buffer;
		for (int pos=0; pos < line_buffer.length(); pos++)
		{
			if (isspace(line_buffer[pos]))
			{
				if (line_buffer[pos] != '\t' && line_buffer[pos] != '\n')
				{
					cerr << "Warning: GFF file has invalid format. A whitespace character other than a tab or newline was detected.\n";
					cerr << "Line producing warning: " << line_buffer << endl;
					break;
				}
			}
		}

		
		// GFF file should have 9 columns separated by tabs
		int i=0;
		int begin = -1;
		int end = -1;
		string id;
		while (getline(stream, field_buffer, '\t'))
		{

			// Fourth column is start of ORF
			if (i == 3)
			{
				begin = stoi(field_buffer);
			}
			// Fifth column is end of ORF
			else if (i == 4)
			{
				end = stoi(field_buffer);
			}

			if (i == 8)
				id = field_buffer;
			
			if (i > 8)
				break;
			
			i++;
		}
		stream.clear();

		// i=0 happens at end of line... do not report error for this
		if (i > 9 && i != 0)
		{
			cerr << "Warning: GFF file has invalid number of columns(" << i+1 << ")\n";
			cerr << "Line producing error: " << line_buffer << endl;
		}

		if (begin != -1 && end != -1)
		{
			
			paramcontainer->cords.push_back(ORF_Cordinate(begin,end, id));
		}
	}

	f.close();
	
	return true;
}

bool get_next_sequences (ifstream *file, vector<string> *queue)
{


	for (int i=0; i < queue->size(); i++)
	{
		string buffer;
	
		do
		{
			if (file->eof())
			{
				// Wipe all the reads in the queue which are old
				for (; i < queue->size(); i++)
					(*queue)[i] = "";
				return true;
			}

			getline(*file, buffer);


			
		} while ((buffer.front() == '>') || !(is_valid_read(buffer)));


		(*queue)[i] = buffer;	
	}

return false;



}

int reverse_read(string *read)
{
	reverse(read->begin(), read->end());
	return 0;
}

int reverse_compliment(string *read)
{
	reverse(read->begin(), read->end());

	for (int i=0; i < read->length(); i++)
	{
		char c = (*read)[i];
		if (c == 'A')
			(*read)[i] = 'T';
		else if (c == 'T')
			(*read)[i] = 'A';
		else if (c == 'G')
			(*read)[i] = 'C';
		else if (c == 'C')
			(*read)[i] = 'G';
		else
		{
			cerr << "Warning: Error creating reverse compliment of read\n";
			return -1;
		}
	}
	return 0;
}


int location_in_genome (string read, string genome)
{
	using namespace seqan;
	int location = -1;

	CharString genome_string = genome;
	CharString read_string = read;

	Finder<String<char>> finder(genome_string);
	Pattern<String<char>, MyersUkkonen> pattern(read_string);

	setScoreLimit(pattern, -1);
	while (find(finder, pattern))
		while (findBegin(finder, pattern, getScore(pattern)))
		{
			int score = getBeginScore(pattern);
			if (score == 0)
			{
				// Substract read length to get real beginning location
				return beginPosition(finder)-read.length();
			}
			else
			{
				// Subtract read length to get real beginning location
				location = beginPosition(finder)-read.length();
			}
		}
	
	return location;	
}


bool is_in_genome (string read, string genome)
{
	if (genome.find(read) != string::npos)
		return true;
	
	return false;
}

string get_genome_chunk (int loc, int read_len, string genome, bool is_circular)
{
	if (loc > genome.length())
	{
		cerr << "Warning: Invalid location given to get_genome_chunk.\n";
		return "";
	}
	
	if (!is_circular)
		return genome.substr(loc, read_len);
	
	if (loc+read_len < genome.length())
		return genome.substr(loc, read_len);
	else
	{		
		// Connect the end and beginning of genome by making this chunk span them				
		int first_part_len = genome.length() - loc;
		int second_part_len = read_len - first_part_len;

		string first = genome.substr(loc, first_part_len);
		string second = genome.substr(0, second_part_len);

		return first+second;
	}
}

/* revise_splice_loc (string read, string genome, int splice_site)
 * Takes a read, genome, and general splice_site as arugments.
 * This function revises the splice_site by finding a match for the
 * first part of the read in the genome.
 */
int revise_splice_site (string read, string genome, int splice_site)
{
	string part1 = read.substr(0, splice_site);
	
	for (int i=0; !(is_in_genome(part1, genome)) && (splice_site-i >= 0); i++)
	{
		part1 = read.substr(0, splice_site-i);
	}

	return part1.length();
}

int array_sum(int array[], int length)
{
    int i;
    int ret = 0;
    for (i=0; i < length; ++i)
        ret = ret + array[i];
    return ret;
}

int determine_splice_loc (string seq, string chunk, struct ParamContainer paramcontainer)
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

    

    for (i=0; i < seq.length(); i++)
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

	
bool splice_is_valid (int splice_site, string part1, string part2, int part1_loc, int part2_loc, string genome, struct ParamContainer paramcontainer)
{
	
	int intron_len = part2_loc - (part1_loc + part1.length());
	
	if (paramcontainer.intron_min > 0 && intron_len < paramcontainer.intron_min)
		return false;
	else if (paramcontainer.intron_max > 0 && intron_len > paramcontainer.intron_max)
		return false;
	

	if (paramcontainer.ignore_noncoding)
	{
		bool found1 = false;
		bool found2 = false;

		for (int i=0; i < paramcontainer.cords.size(); i++)
		{
			int begin, end;
			begin = paramcontainer.cords[i].begin;
			end = paramcontainer.cords[i].end;

			if (part1_loc > begin && part1_loc < end)
			{
				found1 = true;
			}
			
		   	if (part2_loc > begin && part2_loc < end)
			{
				found2 = true;
			}
		}
		
		if (!found1 || !found2)
			return false;
	}


	if (part1.length() > paramcontainer.part1_small_size_threshold)
		if (is_in_genome(part2, genome) && is_in_genome(part1, genome))
				if (part1_loc < part2_loc)
			return true;

	return false;
}

void output_CSV_header ()
{
	cout << "read, donor_loc, acceptor_loc, splice_loc, read_len, intron_len, read_pol, part1_ORF, part2_ORF, same_ORF\n";
}

void output_detected_splice (string read, string part1, int part1_loc, int part2_loc, int splice_loc, int status, ParamContainer paramcontainer)
{
	char read_pol = '+';
	char same_ORF = '-';
	int donor_loc = part1_loc + part1.length();
	int acceptor_loc = part2_loc;
	int intron_len = acceptor_loc - donor_loc;
	int read_len = read.length();
	string part1_ORF = ".";
	string part2_ORF = ".";

	for (int i=0; i < paramcontainer.cords.size(); i++)
	{
		int begin, end;
		begin = paramcontainer.cords[i].begin;
		end = paramcontainer.cords[i].end;

		if (part1_loc > begin && part1_loc < end)
			part1_ORF = paramcontainer.cords[i].id;
		if (part2_loc > begin && part2_loc < end)
			part2_ORF = paramcontainer.cords[i].id;
	}

	if ((part1_ORF != "." && part2_ORF != ".") && (part1_ORF == part2_ORF))
		same_ORF = '+';

	if (status == ON_REVERSE)
	{
		int tmp_loc;

		// Give donor and acceptor locations in relationship to original read
		tmp_loc = donor_loc;
		donor_loc = acceptor_loc;
		acceptor_loc = tmp_loc;

		// Return string to original orientation
		splice_loc = read_len - splice_loc;
	}

	else if (status == ON_ANTI)
	{
		int tmp_loc;

		tmp_loc = donor_loc;
		donor_loc = acceptor_loc;
		acceptor_loc = tmp_loc;

		read_pol = '-';

		splice_loc = read_len - splice_loc;
	}



	cout_mutex.lock();
	cout << read << "," << donor_loc << "," << acceptor_loc << "," << splice_loc << "," << read_len << "," << intron_len << "," << read_pol << "," << part1_ORF << "," << part2_ORF << "," << same_ORF << endl;
	cout_mutex.unlock();
}

void check_for_splice(string read, string genome, struct ParamContainer paramcontainer)
{
	int status = ON_DEFAULT;
	bool finished = false;
	int read_len = read.length();
	string original_read = read;
	
	while (!finished)
	{
		for (int i=0; i < genome.length(); i++)
		{
			if ((!paramcontainer.circular_genome) && (i < genome.length()-read_len))
				break;

			string part1, part2;
			int part1_loc, part2_loc;
			string genome_chunk = get_genome_chunk(i, read_len, genome, paramcontainer.circular_genome);
			int splice_loc = determine_splice_loc(read, genome_chunk, paramcontainer);
			splice_loc = revise_splice_site (read, genome, splice_loc);
			
			if (!(splice_loc > 0))
			{
				finished = true;
				break;
			}
			if (splice_loc >= read.length())
				continue;
			
			part1 = read.substr(0, splice_loc);
			part2 = read.substr(splice_loc);			

		
			// It is known part1 is in the genome_chunk, so only consider it and then add current location(i) to the location
			part1_loc = location_in_genome (part1, genome_chunk);
			part1_loc = part1_loc + i;
		
			part2_loc = location_in_genome (part2, genome);

			if (splice_is_valid (splice_loc, part1, part2, part1_loc, part2_loc, genome, paramcontainer))
			{
				output_detected_splice(original_read, part1, part1_loc, part2_loc, splice_loc, status, paramcontainer);
				finished = true;
				break;
			}


		}
		
		if ((paramcontainer.do_reverse) && (status != ON_REVERSE))
		{
			reverse_read(&read);
			status = ON_REVERSE;
		}
		if ((paramcontainer.do_anti_sense) && (status != ON_ANTI))
		{

			reverse_compliment(&read);
			status = ON_ANTI;
		}
		if ((status == ON_DEFAULT) || (status == ON_ANTI))
			finished = true;
	}
}

int main(int argc, char *argv[])
{
	string filename = "genome.fa";
	string genome;
	string gff_file_loc;
	ifstream read_file;
	bool done = false;
	int threads=0;
	

	struct ParamContainer paramcontainer;
    
    paramcontainer.seed_len = 3;
    paramcontainer.part1_small_size_threshold = 3;
    paramcontainer.part2_small_size_threshold = 3;
    paramcontainer.mM1_ratio_roof = 0.5;
    paramcontainer.mM2_ratio_floor = 0.2;
    paramcontainer.mM1_s_threshold = 3;
    
    paramcontainer.last_n_length = 4;
    paramcontainer.last_n_mM_threshold = 0.20;

	paramcontainer.intron_min = 0;
	paramcontainer.intron_max = 0;

	paramcontainer.do_anti_sense = 1;
	paramcontainer.do_reverse = 1;

	paramcontainer.ignore_noncoding = false;
	paramcontainer.circular_genome = true;

	int queue_size = 10;

	if (argc < 2)
	{
		cerr << "Incorrect number of arguments, consult program documentation\n";
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
		// intron_minimum: --intron_min
		// intron_maximum: --intron_max
		// ignore_noncoding: --ignore-noncoding
		// linear_genome: --linear_genome

		if ( (strncmp(argv[i], "-", 1) == 0) || (strncmp(argv[i], "--", 2) == 0))
		{
			if ( (strcmp(argv[i], "--input") == 0) || (strcmp(argv[i], "-i") == 0))
			{
				cerr << "Input file: " << argv[i+1] << endl;
				read_file.open(argv[i+1]);
				if (!read_file)
				{
					cerr << "Error opening input file located at " << argv[i+1] << endl;
					exit (1);
				}
			}
			else if ( (strcmp(argv[i], "--genome") == 0) || (strcmp(argv[i], "-g") == 0))
			{
				cerr << "Genome file: " << argv[i+1] << endl;
				get_genome(argv[i+1], &genome);
			}
			else if ( (strcmp(argv[i], "--status_log") == 0) || (strcmp(argv[i], "-l") == 0))
			{
				cout << "Status log: " << argv[i+1];
				/*status_log_file_loc = malloc(sizeof(char)*strlen(argv[i+1]));
				strncpy(status_log_file_loc, argv[i+1], strlen(argv[i+1]));*/
			}
			else if ( (strcmp(argv[i], "--threads") == 0) || (strcmp(argv[i], "-t") == 0))
			{
				threads = atoi(argv[i+1]);
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
					cerr << "Warning: --mM1_roof must be between zero and one\n";
			}
			else if (strcmp(argv[i], "--mM2_floor") == 0)
			{
				if ((atof(argv[i+1]) > 0) && (atof(argv[i+1]) < 1))
					paramcontainer.mM2_ratio_floor = atof(argv[i+1]);
				else
					cerr << "Warning: --mM2_floor must be between zero and one\n";
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
			else if (strcmp(argv[i], "--skip-anti-sense") == 0)
			{
				paramcontainer.do_anti_sense = 0;
			}
			else if (strcmp(argv[i], "--skip-reverse") == 0)
			{
				paramcontainer.do_reverse = 0;
			}
			else if (strcmp(argv[i], "--intron_min") == 0)
			{
				paramcontainer.intron_min = atoi(argv[i+1]);
			}
			else if (strcmp(argv[i], "--intron_max") == 0)
			{
				paramcontainer.intron_max = atoi(argv[i+1]);
			}
			else if (strcmp(argv[i], "--ignore-noncoding") == 0)
			{
				paramcontainer.ignore_noncoding = true;
			}
			else if (strcmp(argv[i], "--linear_genome") == 0)
			{
				paramcontainer.circular_genome = false;
			}
			else if (strcmp(argv[i], "-gff") == 0)
			{
				gff_file_loc = argv[i+1];
			}
			else if ((strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "--help") == 0))
			{
				cerr << "Usage: splice_detection -i read_file -g genome_file [options]\n";
				cerr << "See man splice_detection for further information\n";
				return 0;
			}
			else if ((strcmp(argv[i], "-v") == 0) || (strcmp(argv[i], "--version") == 0))
			{
				cerr << "splice_detection " << current_version << endl;
				return 0;
			}
			else
			{
				cerr << "Warning ignoring unknown command: " << argv[i] << endl;
			}
		}

	}


	if (gff_file_loc != "")
	{
		get_ORFS(gff_file_loc, &paramcontainer);
	}
	

	vector<string> read_queue(queue_size);

	
	output_CSV_header();
	
	ThreadPool pool(threads);
	
	do
	{
		done = get_next_sequences(&read_file, &read_queue);

		vector<thread> threads;


		for (int i=0; i < queue_size; i++)
		{
			function<void()> job = bind(check_for_splice, read_queue[i], genome, paramcontainer);
			pool.add_work(job);
		}

	} while (!done);

	pool.terminate();


	return 0;
}
