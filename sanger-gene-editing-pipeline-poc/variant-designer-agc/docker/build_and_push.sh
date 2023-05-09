#!/usr/bin/env bash

# Copyright 2023 Amazon.com and its affiliates; all rights reserved. This file is Amazon Web Services Content and may
# not be duplicated or distributed without permission.

# This script shows how to build the Docker image and push it to ECR to be ready for use
# by Amazon Genomics CLI.

# The argument to this script is the image name. This will be used as the image on the local
# machine and combined with the account and region to form the repository name for ECR.
image=$1
version=$2

if [ "$image" == "" ]
then
    echo "Usage: $0 <image-name> <image-version>"
    exit 1
fi

if [ "$version" == "" ]
then
    echo "Usage: $0 <image-name> <image-version>"
    exit 1
fi


# Get the account number associated with the current IAM credentials
account=$(aws sts get-caller-identity --query Account --output text)

if [ $? -ne 0 ]
then
    exit 255
fi


# Get the region defined in the current configuration (default to us-west-2 if none defined)
region=$(aws configure get region)
region=${region:-us-west-2}

parent_tag="gene-editing-pipeline"
dockerfile="${image}/${version}/Dockerfile"
ecr_repo="${account}.dkr.ecr.${region}.amazonaws.com"
tag_without_version="${parent_tag}/${image}"
tag_with_version="${tag_without_version}:${version}"
local_tag=$tag_with_version
ecr_tag="${ecr_repo}/${tag_with_version}"


echo "AWS Region: ${region}"
echo "Dockerfile: ${dockerfile}"
echo "Local tag : ${local_tag}"
echo "ECR tag   : ${ecr_tag}"


# If the repository doesn't exist in ECR, create it.
aws ecr describe-repositories --repository-names $tag_without_version > /dev/null 2>&1

if [ $? -ne 0 ]
then
    echo "The repository with name ${tag_without_version} does not exist in the registry ${ecr_repo}. Creating repository."
    aws ecr create-repository --repository-name $tag_without_version 
# > /dev/null
fi

# Get the login command from ECR and execute it directly
aws ecr get-login-password --region "${region}" | docker login --username AWS --password-stdin "${account}".dkr.ecr."${region}".amazonaws.com

# Build the docker image locally with the image name and then push it to ECR
# with the full name.

docker buildx build --platform linux/amd64  -t ${local_tag} -f ${dockerfile} ./
docker tag ${local_tag} ${ecr_tag}

docker push ${ecr_tag}

