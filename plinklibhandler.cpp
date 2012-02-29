/*
 * =====================================================================================
 *
 *       Filename:  plinklibhandler.cpp
 *
 *    Description:  Initialization of common plink functions and data structures for use
 *                  in the library.
 *
 *        Created:  02/07/2012
 *
 *         Author:  Nick Davis, nick-davis@utulsa.edu
 *
 * =====================================================================================
 */

#include "options.h"
#include "plink.h"
#include "crandom.h"
#include "perm.h"
#include "helper.h"
#include "sets.h"
#ifdef WIN
#include "time.h"
#endif

#include "plinklibhandler.h"

// global Plink object
Plink* PP;
// log filestream
ofstream LOG;
map<string, int> Range::groupNames;

using namespace std;

// create Plink class and permutation class, seed random, open log filestream
PlinkHandler::PlinkHandler() {
	// Random seed
	if(par::random_seed == 0)
	  CRandom::srand(time(0));
	else
	  CRandom::srand(par::random_seed);

	PP = &ph;

	perm = new Perm(ph);
	PP->pperm = perm;

	// create chromosome definition	
	defineHumanChromosomes();

	// open Plink log
	LOG.open(string(par::output_file_name + ".log").c_str());

	// set fam file name
	par::famfile = par::fileroot + ".fam";
}


// read binary file
void PlinkHandler::readBedFile() {
	par::read_bitfile = true;
	par::bitfilename = par::pedfile = par::fileroot + ".bed";
	par::bitfilename_map = par::mapfile = par::fileroot + ".bim";
	PP->readBinData();
}

// write binary file
void PlinkHandler::writeBedFile() {
	par::write_bitfile = true;
	PP->write_BITFILE();
	if (par::clist) PP->write_covariates();
}

// read plaintext file
void PlinkHandler::readPedFile() {
	par::read_ped = true;
	par::pedfile = par::fileroot + ".ped";
	par::mapfile = par::fileroot + ".map";
	PP->readData();
}

// read numeric file
void PlinkHandler::readNumFile(string numericfile) {
	par::numeric_file = true;
	par::nlist_filename = numericfile;
	if(!PP->readNumericFile()) error("Error: Problem reading the numeric attributes");
}

// read covariate file
void PlinkHandler::readCovFile(string covarfile) {
	par::covar_file = true;
	par::clist = true;
	par::clist_filename = covarfile;
	if(!PP->readCovListFile()) error("Error: Problem reading the covariates");
}

// read pheno file
void PlinkHandler::readPhenoFile(string phenofile) {
	par::pheno_file = true;
	par::pheno_filename = phenofile;
	if(!PP->readPhenoFile())
		error("Error: Problem reading the alternate phenotype file");
}

// prune individuals with missing phenotypes
void PlinkHandler::pruneInd() {
	par::ignore_phenotypes = false;
	removeMissingPhenotypes(*PP);
}

// read extract file
void PlinkHandler::readExtractFile(string extractfile) {
	par::extract_set = true;
	par::extract_file = extractfile;
	PP->extractExcludeSet(false);
}

// read exclude file
void PlinkHandler::readExcludeFile(string excludefile) {
	par::exclude_set = true;
	par::exclude_file = excludefile;
	PP->extractExcludeSet(true);
}

// read remove file
void PlinkHandler::readRemoveFile(string removefile) {
	par::remove_indiv = true;
	par::remove_indiv_list = removefile;
	PP->removeIndividuals(false);
}

// read keep file
void PlinkHandler::readKeepFile(string keepfile) {
	par::keep_indiv = true;
	par::keep_indiv_list = keepfile;
	PP->removeIndividuals(true);
}

// filter founders
void PlinkHandler::filterFounders() {
	par::filter_founders = true;
	PP->filterOnFounder();
}

// initialize Plink data vars, check for dupes, filter SNPs
void PlinkHandler::initData() {
	// Set number of individuals
	PP->n = PP->sample.size();

	// Set number of pairs
	PP->np = (int) ((double) (PP->n * (PP->n - 1)) / (double) 2);

	// Total number of all (test+background) loci
	PP->nl_all = PP->locus.size();

	// Number of permutations to store
	PP->maxr2.resize(par::replicates);

	// Check for duplicate individual or SNP names
	checkDupes(*PP);

	// Binary affection status coding
	if(par::bt)
		affCoding(*PP);

	// Determine formats for printing
	PP->prettyPrintLengths();

	// frequency and genotype filtering
	PP->printLOG("Before frequency and genotyping pruning, there are "
		+ int2str(PP->nl_all) + " SNPs\n");
	PP->filterSNPs();
	if (par::af_write) shutdown();
	PP->printLOG("After frequency and genotyping pruning, there are "
		+ int2str(PP->nl_all) + " SNPs\n");

	// reprint counts
	summaryBasics(*PP);

	// Any null allele codes (monomorhpics)?
	for(int l = 0; l < PP->nl_all; l++) {
		if(PP->locus[l]->allele1 == "")
			PP->locus[l]->allele1 = "0";
	}

	// Set statistics
	S = new Set(PP->snpset);
	PP->pS = S;

	// marker scaffold
	makeScaffold(*PP);
}

// set individual major mode
void PlinkHandler::setInd() {
	PP->SNP2Ind();
}

// perform association test
void PlinkHandler::assocTest() {
	PP->calcAssociationWithPermutation(*PP->pperm);
}

// LD-based pruning
void PlinkHandler::LDPrune() {
	par::prune_ld = true;
	par::prune_ld_pairwise = true;
	par::prune_ld_win = 50;
	par::prune_ld_step = 5;
	par::prune_ld_vif = 0.5;

	PP->pruneLD();
}

// LD statistics
void PlinkHandler::LDStats() {
	PP->calcLDStatistics();
}

// delete encoding data structures
void Plink::cleanUp() {
  for(int i = 0; i < sample.size(); i++)
    delete sample[i];

  for(int l = 0; l < locus.size(); l++)
    delete locus[l];

  if(par::SNP_major)
    for(int l = 0; l < SNP.size(); l++)
      delete SNP[l];
}
