def make_rxn_constraint_string(reaction_ids, lower_bound, upper_bound):
    or_split = [x.strip(' ') for x in reaction_ids.split(' or ')]
    
    if len(or_split) > 1:
        # Handle 'or' case
        reactions = []
        for index, reaction in enumerate(or_split):
            is_reverse = reaction.startswith('reverse_')
            reaction = reaction.replace('reverse_', '')
            if is_reverse:
                reactions.append(f' - {reaction}')
            elif index != 0:
                reactions.append(f' + {reaction}')
            else:
                reactions.append(reaction)
        reactions_str = ''.join(reactions)
        constraint_string = f'{reactions_str} >= {lower_bound}, {reactions_str} <= {upper_bound}'


    else:
        and_split = [x.strip(' ') for x in or_split[0].split(' and ')]
        if len(and_split) > 1:
            # Handle 'and' case
            constraints = []
            for reaction in and_split:
                is_reverse = reaction.startswith('reverse_')
                reaction = reaction.replace('reverse_', '')
                if is_reverse:
                    constraints.append(f'{reaction} >= {-1 * upper_bound}, {reaction} <= {-1 * lower_bound}')
                else:
                    constraints.append(f'{reaction} >= {lower_bound}, {reaction} <= {upper_bound}')
            constraint_string = ', '.join(constraints)
        else:
            # Handle single reaction case
            reaction = and_split[0]
            is_reverse = reaction.startswith('reverse_')
            reaction = reaction.replace('reverse_', '')
            if is_reverse:
                constraint_string = f'{reaction} >= {-1 * upper_bound}, {reaction} <= {-1 * lower_bound}'
            else:
                constraint_string = f'{reaction} >= {lower_bound}, {reaction} <= {upper_bound}'
    
    return constraint_string