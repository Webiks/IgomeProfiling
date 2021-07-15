import boto3 
import botocore.exceptions
import datetime
import re


def is_match_name_instance(names, instance):
    for name in names:
        if re.search(name, instance.tags[0]['Value']):
            return True
    return False


def is_stop_instance(instance, types, names):
    if instance.state['Name'] == 'stopped':
        return False
    if not types and not names:
        return True
    is_type = types and instance.instance_type in types
    is_name = names and is_match_name_instance(names, instance)
    return is_type or is_name


def stop_machines(type_machine, name_machines, logger, region_name = 'us-west-2'):
    logger.info('Stop machines...')
    
    try:
        ec2 = boto3.resource('ec2', region_name=region_name)
    except botocore.exceptions.ClientError as error:
        logger.error(f'{datetime.datetime.now()}: Error {error}  - Can\'t connect to ec2 with boto3 for stop the machines')

    list_instances = ec2.instances.all()

    types = type_machine.split(',') if type_machine else '' 
    names = name_machines.split(',') if name_machines else ''
    for instance in list_instances:
        if is_stop_instance(instance, types, names):
            try:
                instance.stop()
                logger.info(f'Stop the instance - Id: {instance.id}, Type: {instance.instance_type}, Public IPv4: {instance.public_ip_address}, State: {instance.state}')
            except botocore.exceptions.ClientError as error:
                logger.error(f'{datetime.datetime.now()}: Error {error} - can\'t stop the instance {instance.id}')
        else:
            logger.info(f'{datetime.datetime.now()}: Skipping on instance -  {instance.id} and not stop it.')
