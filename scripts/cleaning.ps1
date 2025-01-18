# This script removes rows 1, 3-6 to leave just the genes and column names

Get-ChildItem -Path "data/" -Filter "*.tsv" | ForEach-Object {
    $filePath = $_.FullName
    Write-Host "Processing file: $filePath"
    $lines = Get-Content -Path $filePath
    $filteredLines = $lines | Where-Object {
        $index = $lines.IndexOf($_) + 1
        $index -notin 1, 3, 4, 5, 6
    }
    $filteredLines | Set-Content -Path $filePath
    Write-Host "Completed processing: $filePath"
}
Write-Host "All files have been processed!"