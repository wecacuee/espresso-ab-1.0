#include "set.h"

static pcover fsm_simplify(pset_family F)
{
    pcover D, R;
    D = new_cover(0);
    R = complement(cube1list(F));
    F = espresso(F, D, R);
    free_cover(D);
    free_cover(R);
    return F;
}


void disassemble_fsm(pPLA PLA, int verbose_mode)
{
    int nin, nstates, nout;
    int before, after, present_state, next_state, i, j;
    pcube next_state_mask, present_state_mask, state_mask, p, p1, last;
    pcover go_nowhere, F, tF;

    /* We make the DISGUSTING assumption that the first 'n' outputs have
     *  been created by .symbolic-output, and represent a one-hot encoding
     * of the next state.  'n' is the size of the second-to-last multiple-
     * valued variable (i.e., before the outputs
     */

    if (cube.num_vars - cube.num_binary_vars != 2) {
	fprintf(stderr,
	"use .symbolic and .symbolic-output to specify\n");
	fprintf(stderr,
	"the present state and next state field information\n");
	fatal("disassemble_pla: need two multiple-valued variables\n");
    }

    nin = cube.num_binary_vars;
    nstates = cube.part_size[cube.num_binary_vars];
    nout = cube.part_size[cube.num_vars - 1];
    if (nout < nstates) {
	fprintf(stderr,
	    "use .symbolic and .symbolic-output to specify\n");
	fprintf(stderr,
	    "the present state and next state field information\n");
	fatal("disassemble_pla: # outputs < # states\n");
    }


    present_state = cube.first_part[cube.num_binary_vars];
    present_state_mask = new_cube();
    for(i = 0; i < nstates; i++) {
	set_insert(present_state_mask, i + present_state);
    }

    next_state = cube.first_part[cube.num_binary_vars+1];
    next_state_mask = new_cube();
    for(i = 0; i < nstates; i++) {
	set_insert(next_state_mask, i + next_state);
    }

    state_mask = set_or(new_cube(), next_state_mask, present_state_mask);

    F = new_cover(10);


    /*
     *  check for arcs which go from ANY state to state #i
     */
    for(i = 0; i < nstates; i++) {
	tF = new_cover(10);
	foreach_set(PLA->F, last, p) {
	    if (setp_implies(present_state_mask, p)) { /* from any state ! */
		if (is_in_set(p, next_state + i)) {
		    tF = sf_addset(tF, p);
		}
	    }
	}
	before = tF->count;
	if (before > 0) {
	    tF = fsm_simplify(tF);
	    /* don't allow the next state to disappear ... */
	    foreach_set(tF, last, p) {
		set_insert(p, next_state + i);
	    }
	    after = tF->count;
	    F = sf_append(F, tF);
	    if (verbose_mode) {
		printf("# state EVERY to %d, before=%d after=%d\n",
			i, before, after);
	    }
	}
    }


    /*
     *  some 'arcs' may NOT have a next state -- handle these
     *  we must unravel the present state part
     */
    go_nowhere = new_cover(10);
    foreach_set(PLA->F, last, p) {
	if (setp_disjoint(p, next_state_mask)) { /* no next state !! */
	    go_nowhere = sf_addset(go_nowhere, p);
	}
    }
    before = go_nowhere->count;
    go_nowhere = unravel_range(go_nowhere,
				cube.num_binary_vars, cube.num_binary_vars);
    after = go_nowhere->count;
    F = sf_append(F, go_nowhere);
    if (verbose_mode) {
	printf("# state ANY to NOWHERE, before=%d after=%d\n", before, after);
    }


    /*
     *  minimize cover for all arcs from state #i to state #j
     */
    for(i = 0; i < nstates; i++) {
	for(j = 0; j < nstates; j++) {
	    tF = new_cover(10);
	    foreach_set(PLA->F, last, p) {
		/* not EVERY state */
		if (! setp_implies(present_state_mask, p)) {
		    if (is_in_set(p, present_state + i)) {
			if (is_in_set(p, next_state + j)) {
			    p1 = set_save(p);
			    set_diff(p1, p1, state_mask);
			    set_insert(p1, present_state + i);
			    set_insert(p1, next_state + j);
			    tF = sf_addset(tF, p1);
			    set_free(p1);
			}
		    }
		}
	    }
	    before = tF->count;
	    if (before > 0) {
		tF = fsm_simplify(tF);
		/* don't allow the next state to disappear ... */
		foreach_set(tF, last, p) {
		    set_insert(p, next_state + j);
		}
		after = tF->count;
		F = sf_append(F, tF);
		if (verbose_mode) {
		    printf("# state %d to %d, before=%d after=%d\n",
			    i, j, before, after);
		}
	    }
	}
    }


    free_cube(state_mask);
    free_cube(present_state_mask);
    free_cube(next_state_mask);

    free_cover(PLA->F);
    PLA->F = F;
    free_cover(PLA->D);
    PLA->D = new_cover(0);

    setdown_cube();
    FREE(cube.part_size);
    cube.num_binary_vars = nin;
    cube.num_vars = nin + 3;
    cube.part_size = ALLOC(int, cube.num_vars);
    cube.part_size[cube.num_binary_vars] = nstates;
    cube.part_size[cube.num_binary_vars+1] = nstates;
    cube.part_size[cube.num_binary_vars+2] = nout - nstates;
    cube_setup();

    foreach_set(PLA->F, last, p) {
	kiss_print_cube(stdout, PLA, p, "~1");
    }
}
