/*
 * =====================================================================================
 *
 *       Filename:  plinklibhandler.h
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
#ifndef __PLINKLIBHANDLER_H__
#define __PLINKLIBHANDLER_H__

#include "plink.h"


class PlinkHandler {
	public:
		// create Plink class and permutation class, seed random, open log filestream
		PlinkHandler();
		// read binary file
		void readBedFile();
		// write binary file
		void writeBedFile();
		// read plaintext file
		void readPedFile();
		// read numeric file
		void readNumFile(string numericfile);
		// initialize Plink data vars, check for dupes, filter SNPs
		void initData();
		// set individual major mode
		void setInd();
		// read covariate file
		void readCovFile(string covarfile);
		// read pheno file
		void readPhenoFile(string phenofile);
		// read extract file
		void readExtractFile(string extractfile);
		// read exclude file
		void readExcludeFile(string excludefile);
		// prune individuals with missing phenotypes
		void pruneInd();
		// read remove file
		void readRemoveFile(string removefile);
		// read keep file
		void readKeepFile(string keepfile);
		// filter founders
		void filterFounders();
		// perform association test
		void assocTest();
		// LD-based pruning
		void LDPrune();
		// LD statistics (R and R^2)
		void LDStats();

	private:
		// class Plink object
		Plink ph;
		// Permutation class
		Perm* perm;
		// Set class
		Set* S;

};

#endif
