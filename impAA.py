from operon.components import *
import sys
import os


class Pipeline(ParslPipeline):
	def description(self):
		return "association analysis pipeline using GENESIS and qqman using Parsl"

	def dependencies(self):
		return [
			'numpy',
			'pandas'
			]

	# refer to Parsl documentation -- default Parsl config if user does not provide
	def parsl_configuration(self):
		return {
			'sites':[
				{
					'site': 'Local_Threads',
					'auth': {'channel' : None},
					'execution': {
						'executor': 'threads',
						'provider': None,
						'max_workers': 4
						}
				}
			],
			'globals':{'lazyErrors': True}
		}


	def configuration(self):
		return {
			'R': {
				'path': 'Full path of R exectuable version to be used',
				'lib': 'Full path to directory of R libraries'
			},
			'imputed_files':{
				'path': 'Full path to directory of .dose.vcf.gz and .info.gz files from imputation'
			},
			'dosageConverter':{
				'path': 'Full path of dosageConverter exectuable to be used'
			}
		}


	def arguments(self, parser):
			parser.add_argument('--outdir', help='Path to output directory')
			parser.add_argument('--prefix', type=str, help='Data name or prefix')
			parser.add_argument('--model', type=str, default='logistic', help='logistic or guassian')
			parser.add_argument('--mac', type=int, default=6, help='Minimum minor allele count for filtering out SNPs for qq and manahattan plots. Ex: --mac 6 translates to keeping SNPs with MAC>6 in order to be considered MAC filtered/cleaned')
			parser.add_argument('--maf', type=float, default=0.05, help='Minor allele frequency cut-off for common and rare variants. Ex: --maf 0.05 translates to common>=0.05>rare')
			parser.add_argument('--chunks', type=int, default=100000, help='number of lines to split/chunk variants of vcf for processing')
			parser.add_argument('--tempPath', type=str, default=os.getcwd(), help='full path to exisiting directory to store temp files, i.e. scratch space')
			# same format as chunky

	
	def pipeline(self, pipeline_args, pipeline_config):
		import datetime
		import logging

		'''
		TO DO: add try and except cases
		'''

		def prepFiles(impute_files, chrm, scratch):
			import re
			import gzip
			import itertools
			import logging
			
			logger = logging.getLogger("prepFiles.chr{}".format(str(chrm)))
			logger.info('Prepping chr{} for header and variant extraction...'.format(str(chrm)))
			## for chr{int}.dose.gz files ## OPTIONAL TO DO: RUN IN PARALLEL WITH .info.gz
			doseLine = 0
			doseHeader = open(os.path.join(str(scratch), 'doseHeader_chr{}'.format(str(chrm))), 'wb')
			doseVariants = open(os.path.join(str(scratch), 'doseVars_chr{}'.format(str(chrm))), 'wb')
			
			with gzip.open(os.path.join(impute_files, "chr{}.dose.vcf.gz".format(str(chrm)))) as getDoseHeader:
				for line in enumerate(getDoseHeader):
					if re.search(b'^#', line[1]):
						logger.info('Working...on dosage header extraction')
						doseHeader.write(line[1])
					else:
						logger.info('Dosage file header extraction 100% complete!' + '\n' + 'Now extracting dose variants...')
						doseLine = line[0]
						doseHeader.flush()
						break

			with gzip.open(os.path.join(impute_files, "chr{}.dose.vcf.gz".format(str(chrm)))) as getDoseVars:
				doseVariants.writelines(list(itertools.islice(getDoseVars, doseLine, None)))
				doseVariants.flush()
			
			logger.info('Dosage file header and variant extraction 100% complete')

			## for chr{int}.info.gz files ## OPTIONAL TO DO: RUN IN PARALLEL WITH .dose.gz

			infoLine = 0
			infoHeader = open(os.path.join(str(scratch), 'infoHeader_chr{}'.format(str(chrm))), 'wb')
			infoVariants = open(os.path.join(str(scratch), 'infoVars_chr{}'.format(str(chrm))), 'wb')

			with gzip.open(os.path.join(impute_files, "chr{}.info.gz".format(str(chrm)))) as getInfoHeader:
				for line in enumerate(getInfoHeader):
					if re.search(b'^SNP', line[1]):
						logger.info('Working...on info file header extraction')
						infoHeader.write(line[1])
					else:
						logger.info('Info file header extraction 100% complete!' + '\n' + 'Now extracting info variants...')
						infoLine = line[0]
						infoHeader.flush()
						break

			with gzip.open(os.path.join(impute_files, "chr{}.info.gz".format(str(chrm)))) as getInfoVars:
				infoVariants.writelines(list(itertools.islice(getInfoVars, infoLine, None)))
				infoVariants.flush()

			logger.info('Info file header and variant extraction 100% complete')


		'''
		TO DO: add try and except cases
			   add logging and debugging 
		'''
		def chunkFiles(doseHeader, infoHeader, doseVars, infoVars, chrm, chunks, scratch):
			import re
			import gzip
			import itertools
			import math
			import logging


			logger=logging.getLogger("chunkFiles.chr{}".format(str(chrm)))

			chunkDoseFiles = open(os.path.join(str(scratch), 'chr{}_dose_file_chunks.txt'.format(str(chrm))), 'w')
			chunkInfoFiles = open(os.path.join(str(scratch), 'chr{}_info_file_chunks.txt'.format(str(chrm))), 'w')

			chunkDoseFilesIter = list()
			chunkInfoFilesIter = list()

			with open(os.path.join(pipeline_args['tempPath'], 'doseVars_chr{}'.format(str(chrm)))) as f: lines = sum(1 for line in f)
			# number of files expected to be generated
			totChunksExp = math.ceil(lines/int(pipeline_args['chunks']))
			logger.info("Chunking files on chr{}..total expected file chunks are {}".format(str(chrm), str(totChunksExp)))
	

			# chunking of dose files with header added
			with open(doseVars) as dose:
				#enumerate will label the chunk iterations (i) starting at index 1
				#iter(lambda) now makes iterations of chunks into lists that are enumerated and stored as variable i
				for i, lineChunks in enumerate(iter(lambda:list(itertools.islice(dose, chunks)), []), 1):
					chunkDoseFilesIter.append(os.path.join(str(scratch),"chr{}.dose{}.vcf".format(str(chrm), i)))
					#lineChunks is now a list of with x lines/elements (defined by chunks variable)
					with open(os.path.join(str(scratch),"chr{}.dose{}.vcf".format(str(chrm), i)), 'w') as newDoseFile:
						logger.info("Working on chunk number {} of {} in dose vcf files".format(str(i), str(totChunksExp)))
						#for each list that was generated, it now writes a new file based on the chromosome id and
						#enumerated chunk name
						newDoseFile.write(open(doseHeader, 'r').read())
						newDoseFile.writelines(lineChunks)
						newDoseFile.flush()
				logger.info("COMPLETED {} of {} dose variant files".format(str(totChunksExp), str(totChunksExp)))


			# chunking of info files with header added
			with open(infoVars) as info:
				#enumerate will label the chunk iterations (i) starting at index 1
				#iter(lambda) now makes iterations of chunks into lists that are enumerated and stored as variable i
				for i, lineChunks in enumerate(iter(lambda:list(itertools.islice(info, chunks)), []), 1):
					chunkInfoFilesIter.append(os.path.join(str(scratch),"chr{}.cut{}.info".format(str(chrm), i)))
					#lineChunks is now a list of with x lines/elements (defined by chunks variable)
					with open(os.path.join(str(scratch),"chr{}.cut{}.info".format(str(chrm), i)), 'w') as newInfoFile:
						logger.info("Working on chunk number {} of {} in info files".format(str(i), str(totChunksExp)))

						#for each list that was generated, it now writes a new file based on the chromosome id and
						#enumerated chunk name
						newInfoFile.write(open(infoHeader, 'r').read())
						newInfoFile.writelines(lineChunks)
						newInfoFile.flush()
				logger.info("COMPLETED {} of {} info files".format(str(totChunksExp), str(totChunksExp)))
		


			chunkDoseFiles.write('\n'.join(chunkDoseFilesIter))
			chunkDoseFiles.flush()
			chunkDoseFiles.close()
			chunkInfoFiles.write('\n'.join(chunkInfoFilesIter))	
			chunkInfoFiles.flush()
			chunkInfoFiles.close()	
					
	
		def logicValidation():
			import logging
			print('func: logicValidation() Everything is finished!')
			logging.shutdown()



		## -----------------All code below this line is Software and CodeBlock registration and pipeline logic and flow------------------------
		## ------------------------------------------No functions are defined below this line--------------------------------------------------

		logging.basicConfig(filename=os.path.join(pipeline_args['logs_dir'], 'impAA_run_{}.log'.format(datetime.datetime.now().strftime("%Y-%m-%d_%H:%M"))), level=logging.DEBUG)
		all_chrm_chunks_futures = list()
		'''
		TO DO:  for loop through CodeBlock
				add proper outputs as Data() object
		'''
		for chrm in range(21, 23):
			CodeBlock.register(
				func=prepFiles,
				kwargs={
					'impute_files':pipeline_config['imputed_files']['path'],
					'chrm':str(chrm),
					'scratch':pipeline_args['tempPath']
					},
				outputs=[
						Data(os.path.join(pipeline_args['tempPath'], 'doseHeader_chr{}'.format(str(chrm)))).as_output(),
						Data(os.path.join(pipeline_args['tempPath'], 'doseVars_chr{}'.format(str(chrm)))).as_output(),
						Data(os.path.join(pipeline_args['tempPath'], 'infoHeader_chr{}'.format(str(chrm)))).as_output(),
						Data(os.path.join(pipeline_args['tempPath'], 'infoVars_chr{}'.format(str(chrm)))).as_output()
						]
			)			

		
		'''
		TO DO:  for loop through CodeBlock
		'''

		for chrm in range(21, 23):

			all_chrm_chunks_futures.append(CodeBlock.register(
				func=chunkFiles,
				kwargs={
					'doseHeader':os.path.join(pipeline_args['tempPath'], 'doseHeader_chr{}'.format(str(chrm))), 
					'infoHeader':os.path.join(pipeline_args['tempPath'], 'infoHeader_chr{}'.format(str(chrm))),
					'doseVars':os.path.join(pipeline_args['tempPath'], 'doseVars_chr{}'.format(str(chrm))),
					'infoVars':os.path.join(pipeline_args['tempPath'], 'infoVars_chr{}'.format(str(chrm))),
					'chrm':str(chrm),
					'chunks':pipeline_args['chunks'],
					'scratch':pipeline_args['tempPath']
					},
				inputs=[
					Data(os.path.join(pipeline_args['tempPath'], 'doseHeader_chr{}'.format(str(chrm)))).as_input(),
					Data(os.path.join(pipeline_args['tempPath'], 'doseVars_chr{}'.format(str(chrm)))).as_input(),
					Data(os.path.join(pipeline_args['tempPath'], 'infoHeader_chr{}'.format(str(chrm)))).as_input(),
					Data(os.path.join(pipeline_args['tempPath'], 'infoVars_chr{}'.format(str(chrm)))).as_input()
					],
				outputs=[
					Data(os.path.join(pipeline_args['tempPath'], 'chr{}_dose_file_chunks.txt'.format(str(chrm)))).as_output(),
					Data(os.path.join(pipeline_args['tempPath'], 'chr{}_info_file_chunks.txt'.format(str(chrm)))).as_output()
					]
				)
			)


		CodeBlock.register(
			func=logicValidation,
			wait_on=all_chrm_chunks_futures
			)