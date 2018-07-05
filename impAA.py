from operon.components import *

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
			}
		}


	def arguments(self, parser):
			parser.add_argument('--outdir', help='Path to output directory')
			parser.add_argument('--prefix', type=str, help='Data name or prefix')
			parser.add_argument('--model', type=str, default='logistic', help='logistic or guassian')
			parser.add_argument('--mac', type=int, default=6, help='Minimum minor allele count for filtering out SNPs for qq and manahattan plots. Ex: --mac 6 translates to keeping SNPs with MAC>6 in order to be considered MAC filtered/cleaned')
			parser.add_argument('--maf', type=float, default=0.05, help='Minor allele frequency cut-off for common and rare variants. Ex: --maf 0.05 translates to common>=0.05>rare')
			# same format as chunky

	def pipeline(self, pipeline_args, pipeline_config):
		
		# software initiation
		reformat_impute = Software(name='reformat_impute', path='')	# shell script that reformats .dose.vcf.gz and .info.gz files from imputation	



		for chrm in range(1, 23):
			varname = 'var' + str(chrm)
			reformat_impute.register(
				Parameter(str(varname)),
				Parameter(str(chrm)),
				Redirect(stream=Redirect.BOTH, dest=pipeline_args['outdir']+'reformt_impute_chr{}.log'.format(chrm)),
				extra_outputs=[
					Data('{}.cut{}.info'.format(chrm, split_id)),
					Data('{}.dose{}.vcf'.format(chrm, split_id))
				]
			)
			
