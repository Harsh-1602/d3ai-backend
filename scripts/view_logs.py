#!/usr/bin/env python
"""
Script to view the most recent log entries from the D3AI application.
"""
import os
import sys
import argparse
from pathlib import Path

def tail_file(file_path, num_lines=20):
    """
    Display the last n lines of a file, similar to the Unix 'tail' command.
    
    Args:
        file_path: Path to the log file
        num_lines: Number of lines to display (default: 20)
    """
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            # Read all lines and get the last n
            lines = f.readlines()
            last_lines = lines[-num_lines:] if len(lines) >= num_lines else lines
            
            print(f"\n==== Last {len(last_lines)} lines of {file_path} ====\n")
            for line in last_lines:
                print(line.rstrip())
            print("\n")
            
    except FileNotFoundError:
        print(f"Error: Log file '{file_path}' not found.")
        return 1
    except Exception as e:
        print(f"Error reading log file: {e}")
        return 1
    
    return 0

def main():
    parser = argparse.ArgumentParser(description="View D3AI application logs")
    parser.add_argument(
        "--type", "-t", 
        choices=["all", "error"], 
        default="all",
        help="Type of logs to view (all or error)"
    )
    parser.add_argument(
        "--lines", "-n", 
        type=int, 
        default=20,
        help="Number of lines to display"
    )
    parser.add_argument(
        "--follow", "-f", 
        action="store_true",
        help="Follow the log file (like tail -f)"
    )
    
    args = parser.parse_args()
    
    # Get the project root directory (assuming this script is in the scripts/ folder)
    project_root = Path(__file__).parent.parent
    logs_dir = project_root / "logs"
    
    # Check if logs directory exists
    if not logs_dir.exists():
        print(f"Error: Logs directory '{logs_dir}' not found. Run the application first.")
        return 1
    
    # Determine which log file to read
    log_file = logs_dir / ("error.log" if args.type == "error" else "d3ai.log")
    
    if args.follow:
        try:
            # First print existing content
            tail_file(log_file, args.lines)
            
            # Then follow the file for new content
            print(f"Following {log_file}... (Press Ctrl+C to stop)")
            
            # Get initial file size
            file_size = os.path.getsize(log_file)
            
            while True:
                # Check if file size has changed
                new_size = os.path.getsize(log_file)
                if new_size > file_size:
                    with open(log_file, 'r', encoding='utf-8') as f:
                        # Seek to previous position
                        f.seek(file_size)
                        # Read and print new content
                        new_content = f.read()
                        print(new_content, end='')
                    # Update file size
                    file_size = new_size
                
                # Sleep briefly to reduce CPU usage
                import time
                time.sleep(0.5)
                
        except KeyboardInterrupt:
            print("\nStopped following logs.")
            return 0
    else:
        return tail_file(log_file, args.lines)

if __name__ == "__main__":
    sys.exit(main()) 