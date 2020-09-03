// vim: shiftwidth=4:ts=4:expandtab:cindent
/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2018 Simon Sadedin, ssadedin<at>gmail.com
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
package gngs

/**
 * Wraps a standard OptionAccessor but allows for values to be overridden, enabling
 * classes to be used more generically since OptionAccessor is tied closely to 
 * command line.
 * 
 * @author Simon Sadedin
 */
class CliOptions {
    
    @Delegate
    OptionAccessor opts
    
    HashMap<String, Object> overrides = [:]

    @Override
    Object getProperty(String name) {
        overrides.getOrDefault(name, opts?.getProperty(name))
    }
    
    List<String> arguments() {
        if(overrides.containsKey('arguments'))
            return overrides.arguments
        return opts.arguments()
    }
}
