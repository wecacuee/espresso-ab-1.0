/*
 * Revision Control Information
 *
 * $Source$
 * $Author$
 * $Revision$
 * $Date$
 *
 */
/*
 *  Main driver for espresso
 *
 *  Old style -do xxx, -out xxx, etc. are still supported.
 */

#include "set.h"
#include "main.h"		/* table definitions for options */
#include <unistd.h>

static FILE *last_fp;
static int input_type = FD_type;

void getPLA(int opt, int argc, char **argv, int option, pPLA *PLA, int out_type);
void delete_arg(int *argc, register char **argv, int num);
void init_runtime(void);
void backward_compatibility_hack(int *argc, char **argv, int *option, int *out_type);
void runtime(void);
void usage(void);
bool check_arg(int *argc, register char **argv, register char *s);

int main(int argc, char **argv)
{
    int i, j, first, last, strategy, out_type, option;
    pPLA PLA, PLA1;
    pcover F, Fold, Dold;
    pset last1, p;
    cost_t cost;
    bool error, exact_cover;
    long start;
    extern char *optarg;
    extern int optind;

    start = ptime();

    error = FALSE;
    init_runtime();
#ifdef RANDOM
    srandom(314973);
#endif

    option = 0;			/* default -D: ESPRESSO */
    out_type = F_type;		/* default -o: default is ON-set only */
    debug = 0;			/* default -d: no debugging info */
    verbose_debug = FALSE;	/* default -v: not verbose */
    print_solution = TRUE;	/* default -x: print the solution (!) */
    summary = FALSE;		/* default -s: no summary */
    trace = FALSE;		/* default -t: no trace information */
    strategy = 0;		/* default -S: strategy number */
    first = -1;			/* default -R: select range */
    last = -1;
    remove_essential = TRUE;	/* default -e: */
    force_irredundant = TRUE;
    unwrap_onset = TRUE;
    single_expand = FALSE;
    pos = FALSE;
    recompute_onset = FALSE;
    use_super_gasp = FALSE;
    use_random_order = FALSE;
    kiss = FALSE;
    echo_comments = TRUE;
    echo_unknown_commands = TRUE;
    exact_cover = FALSE;	/* for -qm option, the default */

    backward_compatibility_hack(&argc, argv, &option, &out_type);


    /* parse command line options*/
    while ((i = getopt(argc, argv, "D:S:de:o:r:stv:x")) != EOF) {
	switch(i) {
	    case 'D':		/* -Dcommand invokes a subcommand */
		for(j = 0; option_table[j].name != 0; j++) {
		    if (strcmp(optarg, option_table[j].name) == 0) {
			option = j;
			break;
		    }
		}
		if (option_table[j].name == 0) {
		    fprintf(stderr, "%s: bad subcommand \"%s\"\n",
			argv[0], optarg);
		    exit(1);
		}
		break;

	    case 'o':		/* -ooutput selects and output option */
		for(j = 0; pla_types[j].key != 0; j++) {
		    if (strcmp(optarg, pla_types[j].key+1) == 0) {
			out_type = pla_types[j].value;
			break;
		    }
		}
		if (pla_types[j].key == 0) {
		    fprintf(stderr, "%s: bad output type \"%s\"\n",
			argv[0], optarg);
		    exit(1);
		}
		break;

	    case 'e':		/* -eespresso selects an option for espresso */
		for(j = 0; esp_opt_table[j].name != 0; j++) {
		    if (strcmp(optarg, esp_opt_table[j].name) == 0) {
			*(esp_opt_table[j].variable) = esp_opt_table[j].value;
			break;
		    }
		}
		if (esp_opt_table[j].name == 0) {
		    fprintf(stderr, "%s: bad espresso option \"%s\"\n",
			argv[0], optarg);
		    exit(1);
		}
		break;

	    case 'd':		/* -d turns on (softly) all debug switches */
		debug = debug_table[0].value;
		trace = TRUE;
		summary = TRUE;
		break;

	    case 'v':		/* -vdebug invokes a debug option */
		verbose_debug = TRUE;
		for(j = 0; debug_table[j].name != 0; j++) {
		    if (strcmp(optarg, debug_table[j].name) == 0) {
			debug |= debug_table[j].value;
			break;
		    }
		}
		if (debug_table[j].name == 0) {
		    fprintf(stderr, "%s: bad debug type \"%s\"\n",
			argv[0], optarg);
		    exit(1);
		}
		break;

	    case 't':
		trace = TRUE;
		break;

	    case 's':
		summary = TRUE;
		break;

	    case 'x':		/* -x suppress printing of results */
		print_solution = FALSE;
		break;

	    case 'S':		/* -S sets a strategy for several cmds */
		strategy = atoi(optarg);
		break;

	    case 'r':		/* -r selects range (outputs or vars) */
		if (sscanf(optarg, "%d-%d", &first, &last) < 2) {
		    fprintf(stderr, "%s: bad output range \"%s\"\n",
			argv[0], optarg);
		    exit(1);
		}
		break;

	    default:
		usage();
		exit(1);
	}
    }

    /* provide version information and summaries */
    if (summary || trace) {
	/* echo command line and arguments */
	printf("#");
	for(i = 0; i < argc; i++) {
	    printf(" %s", argv[i]);
	}
	printf("\n");
	printf("# %s\n", VERSION);
    }

    /* the remaining arguments are argv[optind ... argc-1] */
    PLA = PLA1 = NIL(PLA_t);
    switch(option_table[option].num_plas) {
	case 2:
	    if (optind+2 < argc) fatal("trailing arguments on command line");
	    getPLA(optind++, argc, argv, option, &PLA, out_type);
	    getPLA(optind++, argc, argv, option, &PLA1, out_type);
	    break;
	case 1:
	    if (optind+1 < argc) fatal("trailing arguments on command line");
	    getPLA(optind++, argc, argv, option, &PLA, out_type);
	    break;
    }
    if (optind < argc) fatal("trailing arguments on command line");

    if (summary || trace) {
	if (PLA != NIL(PLA_t)) PLA_summary(PLA);
	if (PLA1 != NIL(PLA_t)) PLA_summary(PLA1);
    }

/*
 *  Now a case-statement to decide what to do
 */

    switch(option_table[option].key) {


/******************** Espresso operations ********************/

    case KEY_qm:
	Fold = sf_save(PLA->F);
	PLA->F = minimize_exact(PLA->F, PLA->D, PLA->R, exact_cover);
	EXECUTE(error=verify(PLA->F,Fold,PLA->D), VERIFY_TIME, PLA->F, cost);
	if (error) {
	    print_solution = FALSE;
	    PLA->F = Fold;
	    (void) check_consistency(PLA);
	}
	free_cover(Fold);
	break;

    case KEY_primes:            /* generate all prime implicants */
	EXEC(PLA->F = primes_consensus(cube2list(PLA->F, PLA->D)), 
						    "PRIMES     ", PLA->F);
	break;

    case KEY_map:               /* print out a Karnaugh map of function */
	map(PLA->F);
	print_solution = FALSE;
	break;

/******************** Simple cover operations ********************/

    case KEY_echo:				/* echo the PLA */
	break;

    case KEY_taut:				/* tautology check */
	printf("ON-set is%sa tautology\n",
	    tautology(cube1list(PLA->F)) ? " " : " not ");
	print_solution = FALSE;
	break;

    case KEY_contain:				/* single cube containment */
	PLA->F = sf_contain(PLA->F);
	break;

    case KEY_intersect:				/* cover intersection */
	PLA->F = cv_intersect(PLA->F, PLA1->F);
	break;

    case KEY_union:				/* cover union */
	PLA->F = sf_union(PLA->F, PLA1->F);
	break;

    case KEY_disjoint:				/* make cover disjoint */
	PLA->F = make_disjoint(PLA->F);
	break;

    case KEY_dsharp:				/* cover disjoint-sharp */
	PLA->F = cv_dsharp(PLA->F, PLA1->F);
	break;

    case KEY_sharp:				/* cover sharp */
	PLA->F = cv_sharp(PLA->F, PLA1->F);
	break;

    case KEY_lexsort:				/* lexical sort order */
	PLA->F = lex_sort(PLA->F);
	break;

    case KEY_stats:				/* print info on size */
	if (! summary) PLA_summary(PLA);
	print_solution = FALSE;
	break;

    case KEY_minterms:				/* explode into minterms */
	if (first < 0 || first >= cube.num_vars) {
	    first = 0;
	}
	if (last < 0 || last >= cube.num_vars) {
	    last = cube.num_vars - 1;
	}
	PLA->F = sf_dupl(unravel_range(PLA->F, first, last));
	break;

    }

    /* Print a runtime summary if trace mode enabled */
    if (trace) {
	runtime();
    }

    /* Print total runtime */
    if (summary || trace) {
	print_trace(PLA->F, option_table[option].name, ptime()-start);
    }

    /* Output the solution */
    if (print_solution) {
	EXECUTE(fprint_pla(stdout, PLA, out_type), WRITE_TIME, PLA->F, cost);
    }

    /* Crash and burn if there was a verify error */
    if (error) {
	fatal("cover verification failed");
    }

    /* cleanup all used memory */
    free_PLA(PLA);
    FREE(cube.part_size);
    setdown_cube();             /* free the cube/cdata structure data */
    sf_cleanup();               /* free unused set structures */
    sm_cleanup();               /* sparse matrix cleanup */

    exit(0);
    return 0;
}


void getPLA(int opt, int argc, char **argv, int option, pPLA *PLA, int out_type)
{
    FILE *fp;
    int needs_dcset, needs_offset;
    char *fname;

    if (opt >= argc) {
	fp = stdin;
	fname = "(stdin)";
    } else {
	fname = argv[opt];
	if (strcmp(fname, "-") == 0) {
	    fp = stdin;
	} else if ((fp = fopen(argv[opt], "r")) == NULL) {
	    fprintf(stderr, "%s: Unable to open %s\n", argv[0], fname);
	    exit(1);
	}
    }
    if (option_table[option].key == KEY_echo) {
	needs_dcset = (out_type & D_type) != 0;
	needs_offset = (out_type & R_type) != 0;
    } else {
	needs_dcset = option_table[option].needs_dcset;
	needs_offset = option_table[option].needs_offset;
    }

    if (read_pla(fp, needs_dcset, needs_offset, input_type, PLA) == EOF) {
	fprintf(stderr, "%s: Unable to find PLA on file %s\n", argv[0], fname);
	exit(1);
    }
    (*PLA)->filename = strdup(fname);
    filename = (*PLA)->filename;
/*    (void) fclose(fp);*/
/* hackto support -Dmany */
    last_fp = fp;
}


void runtime(void)
{
    int i;
    long total = 1, temp;

    for(i = 0; i < TIME_COUNT; i++) {
	total += total_time[i];
    }
    for(i = 0; i < TIME_COUNT; i++) {
	if (total_calls[i] != 0) {
	    temp = 100 * total_time[i];
	    printf("# %s\t%2d call(s) for %s (%2ld.%01ld%%)\n",
		total_name[i], total_calls[i], print_time(total_time[i]),
		    temp/total, (10 * (temp%total)) / total);
	}
    }
}


void init_runtime(void)
{
    total_name[READ_TIME] =     "READ       ";
    total_name[WRITE_TIME] =    "WRITE      ";
    total_name[COMPL_TIME] =    "COMPL      ";
    total_name[REDUCE_TIME] =   "REDUCE     ";
    total_name[EXPAND_TIME] =   "EXPAND     ";
    total_name[ESSEN_TIME] =    "ESSEN      ";
    total_name[IRRED_TIME] =    "IRRED      ";
    total_name[GREDUCE_TIME] =  "REDUCE_GASP";
    total_name[GEXPAND_TIME] =  "EXPAND_GASP";
    total_name[GIRRED_TIME] =   "IRRED_GASP ";
    total_name[MV_REDUCE_TIME] ="MV_REDUCE  ";
    total_name[RAISE_IN_TIME] = "RAISE_IN   ";
    total_name[VERIFY_TIME] =   "VERIFY     ";
    total_name[PRIMES_TIME] =   "PRIMES     ";
    total_name[MINCOV_TIME] =   "MINCOV     ";
}


void subcommands(void)
{
    int i, col;
    printf("                ");
    col = 16;
    for(i = 0; option_table[i].name != 0; i++) {
	if ((col + strlen(option_table[i].name) + 1) > 76) {
	    printf(",\n                ");
	    col = 16;
	} else if (i != 0) {
	    printf(", ");
	}
	printf("%s", option_table[i].name);
	col += strlen(option_table[i].name) + 2;
    }
    printf("\n");
}


void usage(void)
{
    printf("%s\n\n", VERSION);
    printf("SYNOPSIS: espresso [options] [file]\n\n");
    printf("  -d        Enable debugging\n");
    printf("  -e[opt]   Select espresso option:\n");
    printf("                fast, ness, nirr, nunwrap, onset, pos, strong,\n");
    printf("                eat, eatdots, kiss, random\n");
    printf("  -o[type]  Select output format:\n");
    printf("                f, fd, fr, fdr, pleasure, eqntott, kiss, cons\n");
    printf("  -rn-m     Select range for subcommands:\n");
    printf("                d1merge: first and last variables (0 ... m-1)\n");
    printf("                minterms: first and last variables (0 ... m-1)\n");
    printf("                opoall: first and last outputs (0 ... m-1)\n");
    printf("  -s        Provide short execution summary\n");
    printf("  -t        Provide longer execution trace\n");
    printf("  -x        Suppress printing of solution\n");
    printf("  -v[type]  Verbose debugging detail (-v '' for all)\n");
    printf("  -D[cmd]   Execute subcommand 'cmd':\n");
    subcommands();
    printf("  -Sn       Select strategy for subcommands:\n");
    printf("                opo: bit2=exact bit1=repeated bit0=skip sparse\n");
    printf("                opoall: 0=minimize, 1=exact\n");
    printf("                pair: 0=algebraic, 1=strongd, 2=espresso, 3=exact\n");
    printf("                pairall: 0=minimize, 1=exact, 2=opo\n");
    printf("                so_espresso: 0=minimize, 1=exact\n");
    printf("                so_both: 0=minimize, 1=exact\n");
}

/*
 *  Hack for backward compatibility (ACK! )
 */

void backward_compatibility_hack(int *argc, char **argv, int *option, int *out_type)
{
    int i, j;

    /* Scan the argument list for something to do (default is ESPRESSO) */
    *option = 0;
    for(i = 1; i < (*argc)-1; i++) {
	if (strcmp(argv[i], "-do") == 0) {
	    for(j = 0; option_table[j].name != 0; j++)
		if (strcmp(argv[i+1], option_table[j].name) == 0) {
		    *option = j;
		    delete_arg(argc, argv, i+1);
		    delete_arg(argc, argv, i);
		    break;
		}
	    if (option_table[j].name == 0) {
		fprintf(stderr,
		 "espresso: bad keyword \"%s\" following -do\n",argv[i+1]);
		exit(1);
	    }
	    break;
	}
    }

    for(i = 1; i < (*argc)-1; i++) {
	if (strcmp(argv[i], "-out") == 0) {
	    for(j = 0; pla_types[j].key != 0; j++)
		if (strcmp(pla_types[j].key+1, argv[i+1]) == 0) {
		    *out_type = pla_types[j].value;
		    delete_arg(argc, argv, i+1);
		    delete_arg(argc, argv, i);
		    break;
		}
	    if (pla_types[j].key == 0) {
		fprintf(stderr,
		   "espresso: bad keyword \"%s\" following -out\n",argv[i+1]);
		exit(1);
	    }
	    break;
	}
    }

    for(i = 1; i < (*argc); i++) {
	if (argv[i][0] == '-') {
	    for(j = 0; esp_opt_table[j].name != 0; j++) {
		if (strcmp(argv[i]+1, esp_opt_table[j].name) == 0) {
		    delete_arg(argc, argv, i);
		    *(esp_opt_table[j].variable) = esp_opt_table[j].value;
		    break;
		}
	    }
	}
    }

    if (check_arg(argc, argv, "-fdr")) input_type = FDR_type;
    if (check_arg(argc, argv, "-fr")) input_type = FR_type;
    if (check_arg(argc, argv, "-f")) input_type = F_type;
}


/* delete_arg -- delete an argument from the argument list */
void delete_arg(int *argc, register char **argv, int num)
{
    register int i;
    (*argc)--;
    for(i = num; i < *argc; i++) {
	argv[i] = argv[i+1];
    }
}


/* check_arg -- scan argv for an argument, and return TRUE if found */
bool check_arg(int *argc, register char **argv, register char *s)
{
    register int i;
    for(i = 1; i < *argc; i++) {
	if (strcmp(argv[i], s) == 0) {
	    delete_arg(argc, argv, i);
	    return TRUE;
	}
    }
    return FALSE;
}
