# Wellcome Sanger Institute Gene Editing Pipeline on AWS
This is a proof of concept to demonstrate how the variant designer workflow of the Sanger gene editing pipeline can 
be implemented in Nextflow DSL and executed on AWS using 
[Amazon Genomics CLI](https://aws.github.io/amazon-genomics-cli/) (ACG). Please see [prerequisites](https://aws.github.io/amazon-genomics-cli/docs/getting-started/prerequisites/) to run AGC.

The sample code; software libraries; command line tools; proofs of concept; templates; or other related technology
(including any of the foregoing that are provided by our personnel) is provided to you as AWS Content under the
AWS Customer Agreement, or the relevant written agreement between you and AWS (whichever applies). You should not
use this AWS Content in your production accounts, or on production or other critical data. You are responsible for
testing, securing, and optimizing the AWS Content, such as sample code, as appropriate for production grade use based
on your specific quality control practices and standards. Deploying AWS Content may incur AWS charges for creating
or using AWS chargeable resources, such as running Amazon EC2 instances or using Amazon S3 storage.

### Prepare workflow execution through AGC ###

#### 1. Download and install AGC 

Follow the [instructions](https://aws.github.io/amazon-genomics-cli/docs/getting-started/installation) on the Amazon
Genomics CLI website.

#### 2. Activate the AWS account for use with AGC

This will deploy the AGC core infrastructure in your AWS account.

```
agc account activate
```

#### 3. Define a username

```
agc configure email you@youremail.com
```

#### 4. Configure additional S3 buckets (optional)

AGC creates an S3 bucket to store logs and outputs and for input caching. If you want to use a separate bucket for 
resources and inputs this needs to be configured in the `agc-project.yaml`:

```
data:
    - location: s3://<your-bucket-name>
      readOnly: true
```

Please note that AGC can write to the bucket provisioned on account activation. Access to any other buckets is read 
only. If you are not using additional S3 buckets delete the data section from `agc-project.yaml`.

#### 5. Provision resources

Navigate to the S3 console to get the name of the S3 bucket created by AGC. Change to the `variant-designer-agc`
directory and run the `provision_resources.sh` script to upload the required input and resource files to the S3 bucket.

```
./provision_resources.sh s3://<agc-bucket-name>
```

#### 6. Build and push Docker container images to ECR

Change into the `docker` directory:

```
cd docker
```

Run `build_and_push.sh` to build the Docker images:

```
./build_and_push.sh autofill-targeton 1.0.0
./build_and_push.sh generate-bedfiles 1.0.0
./build_and_push.sh retrieve-targetons 1.0.0
```

#### 7. Deploy the AGC context

This will deploy the compute environment to execute workflows in your AWS account. Two contexts are defined in 
`agc-project.yaml`: `ondemand` for execution on on-demand EC2 instances and `spot` for execution on spot instances.

To Deploy the `ondemand` context change to the `variant-designer-agc` directory abd run:

```
agc context deploy --context ondemand 
```

#### 8. Edit the `inputs.json` file as required

The `inputs.json` file defines the workflow parameters used by Nextflow to run the workflow, eg:

```
{
    "gene_name" : "CTCF",
    "chromosome" : "chr16",
    "clinvar_ver": "clinvar",
    "targets" : "s3://<agc-bucket-name>/project/GeneEditingPipeline/inputs/VariantDesigner/CTCF_target_regions_manifest.txt",
    "targetons" : "s3://<agc-bucket-name>/project/GeneEditingPipeline/inputs/VariantDesigner/design_choices_51-63.tsv",
    "output_dir":"s3://<agc-bucket-name>/project/GeneEditingPipeline/outputs/VariantDesigner/",
    "validation_utils" : "s3://<agc-bucket-name>/project/GeneEditingPipeline/resources/VariantDesigner/R/validation_utils.R",
    "glue_file" : "s3://<agc-bucket-name>/project/GeneEditingPipeline/resources/VariantDesigner/mane/MANE.GRCh38.v1.0.summary.txt.gz",
    "update_chroms" : "s3://<agc-bucket-name>/project/GeneEditingPipeline/resources/VariantDesigner/other/update_chroms.txt",
    "chrom_file" : "s3://<agc-bucket-name>/project/GeneEditingPipeline/resources/VariantDesigner/assembly/Homo_sapiens_assembly38_chr16.fa",
    "report_dir": "s3://<agc-bucket-name>/project/GeneEditingPipeline/reports/VariantDesigner/",
    "resource_dir": "s3://<agc-bucket-name>/project/GeneEditingPipeline/resources/VariantDesigner",
    "container_registry": "<account-number>.dkr.ecr.us-east-1.amazonaws.com"
}
```

### Execute and track workflows through AGC 

#### 1. Submit a workflow run

```
agc workflow run VariantDesigner -c ondemand
```

#### 2. Check workflow status

```
agc workflow status -c ondemand -r <workflow-instance-id>
```

#### 3. Check Nextflow engine logs

```
agc logs engine -c ondemand -r <workflow-instance-id>
```

#### 4. Check workflow logs

```
agc logs workflow VariantDesigner -r <workflow-instance-id>
```

#### 5. Stop a workflow run 

```
agc workflow stop <workflow-instance-id>
```

See the [AGC command reference](https://aws.github.io/amazon-genomics-cli/docs/reference/) for all agc commands.

### Clean up ### 

#### 1. Destroy the context

This will remove the resources associated with the named context from your account but will keep any S3 outputs and CloudWatch logs.

```
agc context destroy ondemand
```

#### 2. Deactivate the account

If you want stop using Amazon Genomics CLI in your AWS account entirely and remove all resources created by AGC you need to deactivate it.

```
agc account deactivate
```
