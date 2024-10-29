use crate::types::Location;

// merge locations/ranges that abut or overlap
pub fn merge_locations(locations: &mut Vec<Location>)
{
    if locations.len() <= 1 {
        return;
    }

    locations.sort_by_key(|a| a.start);

    for i in (1..locations.len()).rev() {
        eprintln!("{:?}", locations[i]);
        let (this_start, this_end) = {
            let this = locations.get(i).unwrap();
            (this.start, this.end)
        };
        let prev = locations.get_mut(i-1).unwrap();

        if this_start <= prev.end + 1 {
            if this_end > prev.end {
                prev.end = this_end;
            }
            locations.remove(i);
        }
    }
}
